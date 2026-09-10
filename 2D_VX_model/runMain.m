%% runMain: Force-Based Closure Strategy Simulation
%
% Follows force_based_closure_strategy.md (Phases F1/F2/F5) on top of the
% Bell's-law force-based T1/T2 gating with Metropolis-style uphill
% acceptance:
%   - Wound-hole elastic resistance (Kw) decaying over tau_Kw
%   - Purse-string / crawl / margin-contractility all recruited via
%     feedback gated on Kw's decay, not fixed timers
%   - T1 (checkT1transitions.m) / T2 (checkWoundIntercalations.m) gated
%     on real per-edge forces (calibrated f_beta_T1/f_beta_T2), with
%     occasional energy-increasing moves allowed (calibrated
%     T_eff_T1/T_eff_T2)
%
% Saves ALL results to: final_results/
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
param.plotForces = false;

param.woundType = 'circular';
param.woundCircular.center = [param.Lx/2, param.Ly/2];
param.woundCircular.radius = 2.0;

%% 2 - Output directory
saveDir = 'testing';
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

% Reference wound area right after ablation. Uses the same conservation
% measure as the per-step trace (getWoundArea.m) so the two are directly
% comparable; verified to agree with the loop-traced polygon area to 4
% decimal places at t=0, before any pinching can occur.
celldata.A = getCellAreas(celldata, param);
celldata.woundArea0 = getWoundArea(celldata, param);
celldata.Kw_current = param.Kw0;
celldata.lambda_current = 0;
celldata.crawl_current = 0;
celldata.ka_wound_factor_current = 1;
celldata.contractility_wound_current = 1;

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
kaFactorLog = zeros(param.Nsteps,1);
contractilityFactorLog = zeros(param.Nsteps,1);
marginCellLog = zeros(param.Nsteps,1);
woundEdgeLog  = zeros(param.Nsteps,1);
nT1Log        = zeros(param.Nsteps,1);
nSealLog      = zeros(param.Nsteps,1);
nDivLog       = zeros(param.Nsteps,1);

%% Video Setup
videoFilename = fullfile(saveDir, 'wound_simulation.mp4');
vidObj = VideoWriter(videoFilename, 'MPEG-4');
vidObj.FrameRate = 30;
open(vidObj);

fprintf('\n--- Starting simulation: %d timesteps ---\n', param.Nsteps);
totalsimtime = tic;

for tstep = 1:param.Nsteps

    [energymat0] = getTissueEnergyClassical(celldata,param);

    %% Force-based closure state updates (see force_based_closure_strategy.md)
    if celldata.woundArea0 > 0
        celldata.Kw_current = param.Kw0 * exp(-tstep * param.deltat / param.tau_Kw);
        kw_fraction_decayed = 1 - celldata.Kw_current / param.Kw0;
    else
        kw_fraction_decayed = 1;
    end

    target_lambda = param.lambda_purse_string * kw_fraction_decayed;
    celldata.lambda_current = celldata.lambda_current + param.deltat * param.k_recruit * (target_lambda - celldata.lambda_current);

    if param.lambda_purse_string > 0
        celldata.crawl_current = param.crawl_force0 * max(0, 1 - celldata.lambda_current / param.lambda_purse_string);
    else
        celldata.crawl_current = param.crawl_force0;
    end

    target_ka_factor = 1 + (param.ka_wound_factor - 1) * kw_fraction_decayed;
    celldata.ka_wound_factor_current = celldata.ka_wound_factor_current + param.deltat * param.k_recruit * (target_ka_factor - celldata.ka_wound_factor_current);

    target_contractility_factor = 1 + (param.contractility_wound - 1) * kw_fraction_decayed;
    celldata.contractility_wound_current = celldata.contractility_wound_current + param.deltat * param.k_recruit * (target_contractility_factor - celldata.contractility_wound_current);

    kaFactorLog(tstep) = celldata.ka_wound_factor_current;
    contractilityFactorLog(tstep) = celldata.contractility_wound_current;

    %% Forces, position update, geometry refresh
    celldata.f = getVertexForcesClassical(celldata,param,tstep);
    celldata   = updateVertexPositions(celldata,param);
    celldata.A = getCellAreas(celldata,param);
    [celldata.P, celldata.EdgeData] = getCellPerimeters(celldata.nCells,celldata.r,celldata.connec,param,0);

    edgeForceData = getEdgeForces(celldata, param, tstep);

    %% Wound intercalation (T2, force + energy gated)
    if param.enableWoundIntercalations
        [celldata, T1flagVec, intercalation_happened, param] = checkWoundIntercalations(celldata, param, T1flagVec, edgeForceData, tstep);
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
    [celldata, param, seal_happened] = checkWoundSealing(celldata, param, tstep);

    % Wound area by conservation of the periodic box -- no topology
    % assumptions, so it survives the boundary pinching as the wound
    % closes (see getWoundArea.m; loop tracing went NaN from that point
    % on even though the mesh was consistent).
    woundAreaLog(tstep) = getWoundArea(celldata, param);
    marginCellLog(tstep) = numel(param.cellIDtoContract);
    woundEdgeLog(tstep)  = size(celldata.woundEdges,1);
    nT1Log(tstep)        = param.nT1;
    nSealLog(tstep)      = param.nSeals;
    nDivLog(tstep)       = param.nDivisions;

    %% Save a wound/tissue snapshot every timestep for video
    figSnap = figure('visible', 'off');
    plot3Dtissue(celldata.nCells, celldata.r, celldata.connec, param);
    hold on;
    if ~isempty(param.cellIDtoContract)
        p_wound = plot3Dtissue(length(param.cellIDtoContract), celldata.r, celldata.connec(param.cellIDtoContract), param);
        set(p_wound(isgraphics(p_wound)), 'FaceColor', 'magenta', 'FaceAlpha', 0.5);
    end
    title(sprintf('t = %d / %d', tstep, param.Nsteps));
    rectangle('Position',[0 0 param.Lx param.Ly]);
    axis equal;
    hold off;
    
    img = print(figSnap, '-RGBImage');
    writeVideo(vidObj, img);
    
    if mod(tstep, 50) == 0
        saveas(figSnap, fullfile(snapshotDir, sprintf('wound_t%04d.png', tstep)));
    end
    close(figSnap);

    if mod(tstep, 100) == 0
        fprintf('  tstep %4d/%d | area=%6.2f%% | margin cells=%3d | edges=%3d | T1=%3d seals=%2d div=%2d | E=%.3f\n', ...
            tstep, param.Nsteps, 100*woundAreaLog(tstep)/celldata.woundArea0, ...
            marginCellLog(tstep), woundEdgeLog(tstep), param.nT1, param.nSeals, param.nDivisions, energymat(tstep));
    end
end

close(vidObj);
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

%% 6b - Margin cell count / wound edge count vs time (key closure metric)
figMargin = figure('visible','off','Position',[50 50 900 700]);
subplot(2,1,1);
yyaxis left
plot(1:param.Nsteps, marginCellLog, 'LineWidth', 2.0); ylabel('Wound-margin cells');
yyaxis right
plot(1:param.Nsteps, 100*woundAreaLog/celldata.woundArea0, 'LineWidth', 2.0); ylabel('Wound area (% of initial)');
xlabel('Timestep'); title('Wound-margin cell count and wound area vs time'); grid on;
subplot(2,1,2);
plot(1:param.Nsteps, nT1Log, 'LineWidth', 1.5, 'DisplayName','bulk T1'); hold on;
plot(1:param.Nsteps, nSealLog, 'LineWidth', 1.5, 'DisplayName','wound seals');
plot(1:param.Nsteps, nDivLog, 'LineWidth', 1.5, 'DisplayName','cell divisions');
xlabel('Timestep'); ylabel('Cumulative events'); title('Force-gated topological events');
legend('Location','northwest'); grid on;
saveas(figMargin, fullfile(saveDir, 'margin_and_events_vs_time.png'));
close(figMargin);

Tclose = table((1:param.Nsteps)', 100*woundAreaLog/celldata.woundArea0, marginCellLog, woundEdgeLog, nT1Log, nSealLog, nDivLog, ...
    'VariableNames', {'timestep','wound_area_pct','margin_cells','wound_edges','cum_T1','cum_seals','cum_divisions'});
writetable(Tclose, fullfile(saveDir, 'closure_metrics.csv'));

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
fprintf('\n===== SUMMARY =====\n');
fprintf('Cells:            %d -> %d  (divisions=%d)\n', 213, celldata.nCells, param.nDivisions);
fprintf('Wound-margin cells: %d -> %d\n', marginCellLog(1), marginCellLog(end));
fprintf('Wound edges:        %d -> %d\n', woundEdgeLog(1), woundEdgeLog(end));
fprintf('Wound area:         100%% -> %.2f%%  (peak %.2f%%)\n', ...
    100*woundAreaLog(end)/celldata.woundArea0, 100*max(woundAreaLog)/celldata.woundArea0);
fprintf('Events: bulk T1=%d, wound seals=%d, divisions=%d\n', param.nT1, param.nSeals, param.nDivisions);
fprintf('Energy: %.4f -> %.4f\n', energymat(1), energymat(end));
fprintf('Min cell area: %.4f\n', min(celldata.A(celldata.A>0)));
