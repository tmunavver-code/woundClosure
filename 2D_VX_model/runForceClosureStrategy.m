%% runForceClosureStrategy: Force-Based Closure Strategy — 1000-step run
%
% Follows force_based_closure_strategy.md (Phases F1/F2/F5) on top of the
% Bell's-law force-based T1/T2 gating with Metropolis-style uphill
% acceptance:
%   - Kw decaying over tau_Kw as both a recruitment clock and wound-hole
%     elastic resistance when enableWoundResistance is true
%   - Purse-string / crawl / margin-contractility all recruited via
%     feedback gated on Kw's decay, not fixed timers
%   - T1 (checkT1transitions.m) / T2 (checkWoundIntercalations.m) gated
%     on real per-edge forces (calibrated f_beta_T1/f_beta_T2), with
%     occasional energy-increasing moves allowed (calibrated
%     T_eff_T1/T_eff_T2)
%
% Saves ALL results to: results/force_closure_<today's date>/
%   - energy_vs_time.png
%   - ka_contractility_every_50_steps.png + .csv (wound-margin ka/rstiff
%     ramp state sampled every 50 timesteps, plus the step-to-step delta)
%   - wound_area_vs_time.png
%   - initial_state.png / final_state.png
%   - workspace.mat

addpath(genpath('../common'))
addpath(genpath('./functions2D'))

close all;
clear;
clc;

rng('default');

%% 1 - Parameters
param.Lx = 15;
param.Ly = 15;
Case = 'fluid';
SingleCellContractionParameters;

param.Nsteps = 1000;
requestedSteps = str2double(getenv('WOUND_NSTEPS'));
if isfinite(requestedSteps) && requestedSteps >= 1
    param.Nsteps = round(requestedSteps);
end
param.plotForces = false;

param.woundType = 'circular';
param.woundCircular.center = [param.Lx/2, param.Ly/2];
param.woundCircular.radius = 2.0;

% Optional, reproducible local parameter sweep; defaults remain in the parameter file.
configFile = strtrim(getenv('WOUND_CONFIG_JSON'));
if ~isempty(configFile)
    overrides = jsondecode(fileread(configFile));
    overrideNames = fieldnames(overrides);
    for overrideIndex = 1:numel(overrideNames)
        overrideName = overrideNames{overrideIndex};
        if ~isfield(param, overrideName)
            error('Unknown simulation parameter: %s', overrideName);
        end
        param.(overrideName) = overrides.(overrideName);
    end
end


%% 2 - Output directory (unique run folder)
todayStr = datestr(now, 'yyyy-mm-dd_HHMMSS');
if param.enableWoundResistance
    resultPrefix = 'force_closure_gradual_intercalation_with_resistance_';
else
    resultPrefix = 'force_closure_gradual_intercalation_';
end
runTag = strtrim(getenv('WOUND_RUN_TAG'));
if ~isempty(runTag)
    resultPrefix = [resultPrefix runTag '_'];
end
saveDir = fullfile('results', 'TAHAFINAL', [resultPrefix todayStr]);
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end
fprintf('Results will be saved to: %s\n', saveDir);

snapshotDir = fullfile(saveDir, 'snapshots');
if ~exist(snapshotDir, 'dir')
    mkdir(snapshotDir);
end

%% 3 - Generate mesh and create wound
celldata = genDataClassical(param);
[celldata, ~] = createWound(celldata, 'circular', param.woundCircular, param);

[woundEdgeCells, woundEdges] = findWoundEdgeCells(celldata);
param.cellIDtoContract = woundEdgeCells;
celldata.woundEdges = woundEdges;
celldata.cellIDtoContract = woundEdgeCells;

% Reference wound area from the ordered geometric boundary. Area deficit
% is logged separately because division and topology changes can make it a
% poor proxy for the actual opening.
celldata.A = getCellAreas(celldata, param);
initialDiagnostics = getWoundDiagnostics(celldata, param);
if ~initialDiagnostics.loopValid || ~isfinite(initialDiagnostics.loopArea)
    error('Initial wound boundary is not a valid simple loop.');
end
celldata.woundArea0 = initialDiagnostics.loopArea;
celldata.woundPerimeter0 = initialDiagnostics.woundPerimeter;
param.initialWoundBoundaryCells = initialDiagnostics.boundaryCellCount;
celldata.Kw_current = 0;
if isfield(param, 'enforceMonotonicEnergy') && param.enforceMonotonicEnergy
    celldata.lambda_current = param.lambda_purse_string;
    celldata.crawl_current = 0;
    celldata.ka_wound_factor_current = param.ka_wound_factor;
    celldata.contractility_wound_current = param.contractility_wound;
else
    celldata.Kw_current = param.Kw0;
    celldata.lambda_current = 0;
    celldata.crawl_current = 0;
    celldata.ka_wound_factor_current = 1;
    celldata.contractility_wound_current = 1;
end

figInit = figure('visible', 'off');
plot3Dtissue(celldata.nCells, celldata.r, celldata.connec, param);
hold on;
p_wound = plot3Dtissue(length(woundEdgeCells), celldata.r, celldata.connec(woundEdgeCells), param);
set(p_wound(isgraphics(p_wound)), 'FaceColor', 'magenta', 'FaceAlpha', 0.5);
title('Initial State with Wound');
rectangle('Position',[0 0 param.Lx param.Ly]);
hold off;
saveas(figInit, fullfile(saveDir, 'initial_state.png'));
close(figInit);

%% 4 - Run simulation
T1flagVec = zeros(celldata.nMasterVertices,1);
T1relaxstepcountVec = zeros(size(T1flagVec));
T1ForceLog = [];

energymat = zeros(param.Nsteps,1);
timemat   = zeros(param.Nsteps,1);
woundAreaLog = zeros(param.Nsteps,1);
areaDeficitLog = zeros(param.Nsteps,1);
kaFactorLog = zeros(param.Nsteps,1);
contractilityFactorLog = zeros(param.Nsteps,1);
boundaryCellCountLog = zeros(param.Nsteps,1);
boundaryEdgeCountLog = zeros(param.Nsteps,1);
woundPerimeterLog = zeros(param.Nsteps,1);
woundLoopAreaLog = NaN(param.Nsteps,1);
maxVertexDegreeLog = zeros(param.Nsteps,1);
maxWoundVertexDegreeLog = zeros(param.Nsteps,1);
minCellAreaLog = NaN(param.Nsteps,1);
minEdgeLengthLog = NaN(param.Nsteps,1);
topologyValidLog = false(param.Nsteps,1);
motionStepScaleLog = ones(param.Nsteps,1);
intercalationEventLog = cell(param.Nsteps,1);
sealEventLog = cell(param.Nsteps,1);

%% Optional every-timestep video recording
% Enable with WOUND_SAVE_VIDEO=1.  The first frame is the initial state
% (timestep 0), followed by one frame after each completed timestep.
saveVideo = any(strcmpi(strtrim(getenv('WOUND_SAVE_VIDEO')), ...
    {'1', 'true', 'yes'}));
videoFilename = '';
figVideo = [];
if saveVideo
    videoFrameRate = str2double(getenv('WOUND_VIDEO_FPS'));
    if ~isfinite(videoFrameRate) || videoFrameRate <= 0
        videoFrameRate = 20;
    end
    videoFilename = fullfile(saveDir, ...
        'wound_closure_every_timestep.mp4');
    vidObj = VideoWriter(videoFilename, 'MPEG-4');
    vidObj.FrameRate = videoFrameRate;
    vidObj.Quality = 95;
    open(vidObj);

    figVideo = figure('visible', 'off', 'Color', 'white', ...
        'Position', [100 100 960 720]);
    plot3Dtissue(celldata.nCells, celldata.r, celldata.connec, param);
    hold on;
    if ~isempty(param.cellIDtoContract)
        p_wound = plot3Dtissue(length(param.cellIDtoContract), ...
            celldata.r, celldata.connec(param.cellIDtoContract), param);
        set(p_wound(isgraphics(p_wound)), 'FaceColor', 'magenta', ...
            'FaceAlpha', 0.5);
    end
    title(sprintf(['Timestep 0 / %d   |   wound area 100.00%%   |   ' ...
        'boundary cells %d'], param.Nsteps, ...
        initialDiagnostics.boundaryCellCount));
    rectangle('Position', [0 0 param.Lx param.Ly]);
    axis equal;
    xlim([0 param.Lx]);
    ylim([0 param.Ly]);
    hold off;
    drawnow;
    writeVideo(vidObj, print(figVideo, '-RGBImage', '-r120'));
    fprintf('Recording every timestep to: %s (%.1f fps)\n', ...
        videoFilename, videoFrameRate);
end

fprintf('\n--- Starting simulation: %d timesteps ---\n', param.Nsteps);
totalsimtime = tic;

for tstep = 1:param.Nsteps

    [energymat0] = getTissueEnergyClassical(celldata,param);

    %% Closure-state update
    if isfield(param, 'enforceMonotonicEnergy') && param.enforceMonotonicEnergy
        % Keep the energy functional fixed in time. There is no wound-area
        % resistance or non-conservative crawling in this mode.
        celldata.Kw_current = 0;
        celldata.lambda_current = param.lambda_purse_string;
        celldata.crawl_current = 0;
        celldata.ka_wound_factor_current = param.ka_wound_factor;
        celldata.contractility_wound_current = param.contractility_wound;
    else
        if celldata.woundArea0 > 0
            celldata.Kw_current = param.Kw0 * exp( ...
                -tstep * param.deltat / param.tau_Kw);
            kw_fraction_decayed = 1 - celldata.Kw_current / param.Kw0;
        else
            kw_fraction_decayed = 1;
        end

        target_lambda = param.lambda_purse_string * kw_fraction_decayed;
        celldata.lambda_current = celldata.lambda_current + ...
            param.deltat * param.k_recruit * ...
            (target_lambda - celldata.lambda_current);

        if param.lambda_purse_string > 0
            celldata.crawl_current = param.crawl_force0 * max(0, ...
                1 - celldata.lambda_current / param.lambda_purse_string);
        else
            celldata.crawl_current = param.crawl_force0;
        end

        target_ka_factor = 1 + ...
            (param.ka_wound_factor - 1) * kw_fraction_decayed;
        celldata.ka_wound_factor_current = ...
            celldata.ka_wound_factor_current + ...
            param.deltat * param.k_recruit * ...
            (target_ka_factor - celldata.ka_wound_factor_current);

        target_contractility_factor = 1 + ...
            (param.contractility_wound - 1) * kw_fraction_decayed;
        celldata.contractility_wound_current = ...
            celldata.contractility_wound_current + ...
            param.deltat * param.k_recruit * ...
            (target_contractility_factor - ...
            celldata.contractility_wound_current);
    end

    kaFactorLog(tstep) = celldata.ka_wound_factor_current;
    contractilityFactorLog(tstep) = celldata.contractility_wound_current;

    %% Forces, position update, geometry refresh
    celldata.f = getVertexForcesClassical(celldata,param,tstep);
    if isfield(param, 'enforceMonotonicEnergy') && param.enforceMonotonicEnergy
        [celldata, motionStepScaleLog(tstep)] = ...
            updateVertexPositionsMonotonic(celldata, param);
    else
        celldata = updateVertexPositions(celldata,param);
        celldata.A = getCellAreas(celldata,param);
        [celldata.P, celldata.EdgeData] = getCellPerimeters( ...
            celldata.nCells, celldata.r, celldata.connec, param, 0);
    end

    edgeForceData = getEdgeForces(celldata, param, tstep);

    %% Wound intercalation (T2, force + energy gated)
    if param.enableWoundIntercalations
        [celldata, T1flagVec, intercalation_happened, param, intercalationEvent] = ...
            checkWoundIntercalations(celldata, param, T1flagVec, edgeForceData, tstep);
        intercalationEventLog{tstep} = intercalationEvent;
        if intercalation_happened
            celldata.A = getCellAreas(celldata,param);
            [celldata.P, celldata.EdgeData] = getCellPerimeters(celldata.nCells,celldata.r,celldata.connec,param,0);
            edgeForceData = getEdgeForces(celldata, param, tstep);
        end
    end

    %% T1 (force + energy gated, Metropolis uphill acceptance)
    nT1_before = param.nT1;
    if param.enableT1transitions
        [celldata.r,celldata.connec,celldata.EdgeData,celldata.verttocell, T1flagVec,T1relaxstepcountVec,param.nT1, T1ForceLog] = ...
            checkT1transitions(celldata.nCells,celldata.r,celldata.connec,celldata.EdgeData,celldata.verttocell,param,T1flagVec,T1relaxstepcountVec,0, edgeForceData, tstep, T1ForceLog, celldata);
    end
    if param.nT1 > nT1_before
        % Bulk T1 rewired connectivity -- the wound margin may have
        % changed. See refreshWoundMargin.m for why staleness here is a
        % real physics error, not just bookkeeping.
        [celldata, param] = refreshWoundMargin(celldata, param);
    end

    param.Tsim = param.Tsim + param.deltat;
    celldata.r0 = celldata.r;

    %% Cell division (force-based, third topological event type)
    % See checkCellDivision.m / performCellDivision.m. Grows
    % nMasterVertices/nCells -- every per-vertex array sized at the old
    % nMasterVertices must grow with it.
    nMasterVerticesBefore = celldata.nMasterVertices;
    [celldata, param, division_happened] = checkCellDivision(celldata, param);
    if division_happened
        nNewVerts = celldata.nMasterVertices - nMasterVerticesBefore;
        if nNewVerts > 0
            T1flagVec(end+1:end+nNewVerts) = 0;
            T1relaxstepcountVec(end+1:end+nNewVerts) = 0;
            celldata.f(end+1:end+nNewVerts, :) = 0;
            celldata.r0(end+1:end+nNewVerts, :) = celldata.r(end-nNewVerts+1:end, :);
        end
        celldata.A = getCellAreas(celldata,param);
        [celldata.P, celldata.EdgeData] = getCellPerimeters(celldata.nCells,celldata.r,celldata.connec,param,0);
        % Division rewired connectivity and added cells -- refresh margin.
        [celldata, param] = refreshWoundMargin(celldata, param);
    end

    energymat(tstep) = getTissueEnergyClassical(celldata,param);
    timemat(tstep)   = param.Tsim;

    %% Wound sealing (force + energy gated) -- retires cells from the
    % margin as opposing sides come into contact. See checkWoundSealing.m.
    [celldata, param, seal_happened, sealEvent] = ...
        checkWoundSealing(celldata, param, tstep);
    sealEventLog{tstep} = sealEvent;
    if seal_happened
        celldata.A = getCellAreas(celldata,param);
        [celldata.P, celldata.EdgeData] = getCellPerimeters( ...
            celldata.nCells, celldata.r, celldata.connec, param, 0);
    end
    energymat(tstep) = getTissueEnergyClassical(celldata,param);

    %% Wound geometry and topology diagnostics
    diagnostics = getWoundDiagnostics(celldata, param);
    woundAreaLog(tstep) = diagnostics.loopArea;
    areaDeficitLog(tstep) = diagnostics.areaDeficitRaw;
    boundaryCellCountLog(tstep) = diagnostics.boundaryCellCount;
    boundaryEdgeCountLog(tstep) = diagnostics.boundaryEdgeCount;
    woundPerimeterLog(tstep) = diagnostics.woundPerimeter;
    woundLoopAreaLog(tstep) = diagnostics.loopArea;
    maxVertexDegreeLog(tstep) = diagnostics.maxVertexDegree;
    maxWoundVertexDegreeLog(tstep) = diagnostics.maxWoundVertexDegree;
    minCellAreaLog(tstep) = diagnostics.minPositiveCellArea;
    minEdgeLengthLog(tstep) = diagnostics.minEdgeLength;
    topologyValidLog(tstep) = diagnostics.topologyValid && diagnostics.loopValid;
    if ~topologyValidLog(tstep)
        error('Invalid topology at step %d: %s (loop valid: %d)', ...
            tstep, diagnostics.topologyReason, diagnostics.loopValid);
    end

    %% Render the state for video and/or the 50-step PNG snapshots
    if saveVideo || mod(tstep, 50) == 0
        if saveVideo
            clf(figVideo);
            figure(figVideo);
            figSnap = figVideo;
        else
            figSnap = figure('visible', 'off', 'Color', 'white', ...
                'Position', [100 100 960 720]);
        end
        plot3Dtissue(celldata.nCells, celldata.r, celldata.connec, param);
        hold on;
        if ~isempty(param.cellIDtoContract)
            p_wound = plot3Dtissue(length(param.cellIDtoContract), celldata.r, celldata.connec(param.cellIDtoContract), param);
            set(p_wound(isgraphics(p_wound)), 'FaceColor', 'magenta', 'FaceAlpha', 0.5);
        end
        title(sprintf(['Timestep %d / %d   |   wound area %.2f%%   |   ' ...
            'boundary cells %d   |   energy %.3f'], ...
            tstep, param.Nsteps, ...
            100*woundAreaLog(tstep)/celldata.woundArea0, ...
            boundaryCellCountLog(tstep), energymat(tstep)));
        rectangle('Position',[0 0 param.Lx param.Ly]);
        axis equal;
        xlim([0 param.Lx]);
        ylim([0 param.Ly]);
        hold off;
        drawnow;
        if saveVideo
            writeVideo(vidObj, print(figSnap, '-RGBImage', '-r120'));
        end
        if mod(tstep, 50) == 0
            saveas(figSnap, fullfile(snapshotDir, ...
                sprintf('wound_t%04d.png', tstep)));
        end
        if ~saveVideo
            close(figSnap);
        end
    end

    if mod(tstep, 100) == 0
        fprintf(['  tstep %d / %d, energy=%.4f, nT1=%d, wound area=%.2f%%, ' ...
            'boundary cells=%d, max degree=%d\n'], tstep, param.Nsteps, ...
            energymat(tstep), param.nT1, ...
            100*woundAreaLog(tstep)/celldata.woundArea0, ...
            boundaryCellCountLog(tstep), maxVertexDegreeLog(tstep));
    end
end

if saveVideo
    close(vidObj);
    close(figVideo);
    fprintf('Video saved with %d frames: %s\n', ...
        param.Nsteps + 1, videoFilename);
end

fprintf('Simulation finished in %.1f s\n', toc(totalsimtime));

%% 5 - Save energy vs time
figEnergy = figure('visible', 'off');
plot(timemat, energymat, 'LineWidth', 2.0);
xlabel('Simulation time'); ylabel('Total tissue energy');
title('Energy vs Time'); grid on;
saveas(figEnergy, fullfile(saveDir, 'energy_vs_time.png'));
close(figEnergy);

%% 6 - Save wound area vs time
figArea = figure('visible', 'off');
plot(1:param.Nsteps, 100*woundAreaLog/celldata.woundArea0, 'LineWidth', 2.0);
xlabel('Timestep'); ylabel('Wound area (% of initial)');
title('Wound Area vs Time'); grid on;
saveas(figArea, fullfile(saveDir, 'wound_area_vs_time.png'));
close(figArea);

%% 6b - Save wound-boundary topology diagnostics
figBoundary = figure('visible', 'off', 'Position', [50 50 900 750]);
subplot(3,1,1);
plot(1:param.Nsteps, boundaryCellCountLog, 'LineWidth', 2.0);
xlabel('Timestep'); ylabel('Boundary cells'); grid on;
title('Cells touching the wound boundary');
subplot(3,1,2);
plot(100*woundAreaLog/celldata.woundArea0, boundaryCellCountLog, ...
    'LineWidth', 2.0);
xlabel('Wound area (% of initial)'); ylabel('Boundary cells'); grid on;
set(gca, 'XDir', 'reverse');
title('Boundary-cell loss during closure');
subplot(3,1,3);
plot(1:param.Nsteps, maxVertexDegreeLog, 'LineWidth', 1.5, ...
    'DisplayName', 'global');
hold on;
plot(1:param.Nsteps, maxWoundVertexDegreeLog, 'LineWidth', 1.5, ...
    'DisplayName', 'wound');
yline(param.maxSealVertexDegree, '--', 'seal limit');
xlabel('Timestep'); ylabel('Maximum incident cells'); grid on;
legend('Location', 'best');
title('Vertex coordination');
saveas(figBoundary, fullfile(saveDir, 'wound_topology_diagnostics.png'));
close(figBoundary);

diagnosticTable = table((1:param.Nsteps)', timemat, woundAreaLog, ...
    areaDeficitLog, ...
    woundPerimeterLog, boundaryCellCountLog, ...
    boundaryEdgeCountLog, maxVertexDegreeLog, maxWoundVertexDegreeLog, ...
    minCellAreaLog, minEdgeLengthLog, topologyValidLog, ...
    'VariableNames', {'timestep','simulation_time','loop_area', ...
    'area_deficit','wound_perimeter','boundary_cells','boundary_edges', ...
    'max_vertex_degree','max_wound_vertex_degree','min_cell_area', ...
    'min_edge_length','topology_valid'});
writetable(diagnosticTable, fullfile(saveDir, 'wound_topology_diagnostics.csv'));

intercalationEvents = vertcat(intercalationEventLog{:});
writetable(struct2table(intercalationEvents), ...
    fullfile(saveDir, 'wound_intercalation_events.csv'));
sealEvents = vertcat(sealEventLog{:});
writetable(struct2table(sealEvents), fullfile(saveDir, 'wound_seal_events.csv'));

%% 7 - Save ka / contractility ramp, sampled every 50 timesteps, with the
%      step-to-step (50-step-block) difference
sampleSteps = 50:50:param.Nsteps;
kaSamples = kaFactorLog(sampleSteps);
contrSamples = contractilityFactorLog(sampleSteps);
kaDiff = [kaSamples(1) - 1; diff(kaSamples)];         % vs. baseline (1.0) for the first sample
contrDiff = [contrSamples(1) - 1; diff(contrSamples)];

figKa = figure('visible', 'off', 'Position', [50 50 900 700]);
subplot(2,1,1);
plot(sampleSteps, kaSamples, '-o', 'LineWidth', 1.5, 'DisplayName', 'ka factor');
hold on;
plot(sampleSteps, contrSamples, '-s', 'LineWidth', 1.5, 'DisplayName', 'contractility factor');
xlabel('Timestep'); ylabel('Wound-margin multiplier (x bulk)');
title('Wound-margin contractility ramp, sampled every 50 steps');
legend('Location','southeast'); grid on;

subplot(2,1,2);
bar(sampleSteps, kaDiff, 'DisplayName', 'ka factor \Delta per 50 steps');
hold on;
plot(sampleSteps, contrDiff, '-r', 'LineWidth', 1.5, 'DisplayName', 'contractility factor \Delta per 50 steps');
xlabel('Timestep'); ylabel('Change since previous 50-step sample');
title('Timestep-50 difference in ka / contractility ramp');
legend('Location','northeast'); grid on;
saveas(figKa, fullfile(saveDir, 'ka_contractility_every_50_steps.png'));
close(figKa);

% CSV of the raw 50-step samples
T = table(sampleSteps', kaSamples, kaDiff, contrSamples, contrDiff, ...
    'VariableNames', {'timestep','ka_factor','ka_factor_delta','contractility_factor','contractility_factor_delta'});
writetable(T, fullfile(saveDir, 'ka_contractility_every_50_steps.csv'));

%% 8 - Final state snapshot
figFinal = figure('visible', 'off');
plot3Dtissue(celldata.nCells, celldata.r, celldata.connec, param);
hold on;
if ~isempty(param.cellIDtoContract)
    p_wound = plot3Dtissue(length(param.cellIDtoContract), celldata.r, celldata.connec(param.cellIDtoContract), param);
    set(p_wound(isgraphics(p_wound)), 'FaceColor', 'magenta', 'FaceAlpha', 0.5);
end
title(sprintf('Final State (t = %d steps)', param.Nsteps));
rectangle('Position',[0 0 param.Lx param.Ly]);
hold off;
saveas(figFinal, fullfile(saveDir, 'final_state.png'));
close(figFinal);

%% 9 - Save workspace
save(fullfile(saveDir, 'workspace.mat'));

fprintf('\nAll results saved to: %s\n', saveDir);
fprintf('Total T1 transitions: %d\n', param.nT1);
fprintf('Wound handoffs: %d; corner trims: %d; seals: %d\n', ...
    param.nWoundIntercalations, param.nWoundCornerTrims, param.nSeals);
fprintf('Final wound area: %.2f%% of initial\n', 100*woundAreaLog(end)/celldata.woundArea0);
fprintf('Maximum vertex degree: %d; all topology checks valid: %d\n', ...
    max(maxVertexDegreeLog), all(topologyValidLog));
