function metrics = runEnergyParameterTrial(configPath,tag)
setenv('WOUND_CONFIG_JSON',configPath);
setenv('WOUND_RUN_TAG',tag);
runClock=tic;
addpath(genpath('../common'))
addpath(genpath('./functions2D'))





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

energyPhases=zeros(param.Nsteps,6);
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

    energyPhases(tstep,1:2) = [energymat0, getTissueEnergyClassical(celldata,param)];
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

    energyPhases(tstep,3) = getTissueEnergyClassical(celldata,param);
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

    energyPhases(tstep,4) = getTissueEnergyClassical(celldata,param);
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

    energyPhases(tstep,5) = getTissueEnergyClassical(celldata,param);
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

    energyPhases(tstep,6) = energymat(tstep);
end

elapsedSeconds=toc(runClock);
energyPhasesTable=array2table([(1:param.Nsteps)',energyPhases], 'VariableNames', ...
 {'step','before_recruitment','after_recruitment','after_motion','after_handoff','after_T1','final'});
writetable(energyPhasesTable,fullfile(saveDir,'energy_budget.csv'));
intercalationEvents=vertcat(intercalationEventLog{:});
writetable(struct2table(intercalationEvents),fullfile(saveDir,'wound_intercalation_events.csv'));
metrics=struct('folder',saveDir,'runtime_seconds',elapsedSeconds,...
 'energy_start',energymat(1),'energy_end',energymat(end),...
 'up_steps',nnz(diff(energymat)>1e-10),'total_up',sum(max(diff(energymat),0)),...
 'max_up',max([0;diff(energymat)]),'area_final_pct',100*woundAreaLog(end)/celldata.woundArea0,...
 'boundary_final',boundaryCellCountLog(end),'handoffs',param.nWoundIntercalations,...
 'T1',param.nT1,'valid',all(topologyValidLog),'max_degree',max(maxVertexDegreeLog));
fid=fopen(fullfile(saveDir,'metrics.json'),'w'); fprintf(fid,'%s',jsonencode(metrics)); fclose(fid);
save(fullfile(saveDir,'workspace.mat'),'-v7');
fid=fopen(fullfile(saveDir,'parameters.json'),'w'); fprintf(fid,'%s',jsonencode(param)); fclose(fid);
f=figure('Visible','off'); plot(timemat,energymat,'LineWidth',1.5); grid on;
xlabel('Simulation time'); ylabel('Total mechanical potential energy');
title(strrep(getenv('WOUND_RUN_TAG'),'_',' ')); saveas(f,fullfile(saveDir,'energy_vs_time.png')); close(f);
f=figure('Visible','off'); plot(timemat,100*woundAreaLog/celldata.woundArea0,'LineWidth',1.5);grid on;
xlabel('Simulation time');ylabel('Wound area (% initial)');saveas(f,fullfile(saveDir,'wound_area_vs_time.png'));close(f);
f=figure('Visible','off');plot(timemat,boundaryCellCountLog,'LineWidth',1.5);grid on;
xlabel('Simulation time');ylabel('Wound boundary cells');saveas(f,fullfile(saveDir,'boundary_cells_vs_time.png'));close(f);
f=figure('Visible','off');plot3Dtissue(celldata.nCells,celldata.r,celldata.connec,param);
axis equal;saveas(f,fullfile(saveDir,'final_state.png'));close(f);
fprintf('SEARCH_RESULT %s\n',jsonencode(metrics));
end
