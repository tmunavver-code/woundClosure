%% runEdgeForceAnalysis: Standalone script to compute and visualize edge forces
%
% This script runs the full 2D vertex model simulation and, at the end,
% computes the force distribution on every cell edge. It produces:
%   1. A 4-panel figure showing perimeter tension, purse-string tension,
%      total tension, and area pressure on each edge.
%   2. A vector plot showing tangential + normal force arrows at edge midpoints.
%   3. Saved results (workspace + figures) to results/edge_forces/
%
% Usage: Simply run this script in MATLAB from the 2D_VX_model directory.

addpath(genpath('../common'))
addpath(genpath('./functions2D'))

close all;
clear;
clc;

rng('default');

%% 1 - Initialising the parameters
param.Lx = 15;
param.Ly = 15;

Case = 'fluid';

SingleCellContractionParameters;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% --- WOUND CONFIGURATION --- %%%
param.woundType = 'circular';
param.woundCircular.center = [param.Lx/2, param.Ly/2];
param.woundCircular.radius = 2.0;
param.woundElliptical.center = [param.Lx/2, param.Ly/2];
param.woundElliptical.x_axis = 8.0;
param.woundElliptical.y_axis = 1.0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 2 - Generating mesh and Creating Wound
celldata = genDataClassical(param);

switch param.woundType
    case 'circular'
        [celldata, ~] = createWound(celldata, 'circular', param.woundCircular, param);
    case 'elliptical'
        [celldata, ~] = createWound(celldata, 'elliptical', param.woundElliptical, param);
    case 'none'
        disp('No wound will be created.');
end

[woundEdgeCells, woundEdges] = findWoundEdgeCells(celldata);
param.cellIDtoContract = woundEdgeCells;
celldata.woundEdges = woundEdges;

%% 3 - Run simulation
totalsimtime = tic;

energymat = zeros(param.Nsteps, 1);
timemat   = zeros(param.Nsteps, 1);
T1flagVec = zeros(celldata.nMasterVertices, 1);
T1relaxstepcountVec = zeros(size(T1flagVec));

Coordinates = zeros(celldata.nMasterVertices, 2*param.Nsteps);
Connectivity = cell(celldata.nCells, param.Nsteps);

% Store edge force data at selected timesteps
edgeForceSnapshots = struct();
snapshotTimesteps = [1, round(param.Nsteps/4), round(param.Nsteps/2), ...
                     round(3*param.Nsteps/4), param.Nsteps];
snapshotCount = 0;

for tstep = 1:param.Nsteps
    
    applyStretchX;
    
    [energymat0] = getTissueEnergyClassical(celldata, param);
    
    celldata.f = getVertexForcesClassical(celldata, param, tstep);
    
    if param.enableWoundIntercalations
        [celldata, T1flagVec, intercalation_happened, param] = checkWoundIntercalations(celldata, param, T1flagVec);
        if intercalation_happened
            celldata.A = getCellAreas(celldata, param);
            [celldata.P, celldata.EdgeData] = getCellPerimeters(celldata.nCells, celldata.r, celldata.connec, param, 0);
        end
    end
    
    celldata = updateVertexPositions(celldata, param);
    
    celldata.A = getCellAreas(celldata, param);
    [celldata.P, celldata.EdgeData] = getCellPerimeters(celldata.nCells, celldata.r, celldata.connec, param, 0);
    
    if param.enableT1transitions
        [celldata.r, celldata.connec, celldata.EdgeData, celldata.verttocell, ...
         T1flagVec, T1relaxstepcountVec, param.nT1] = ...
            checkT1transitions(celldata.nCells, celldata.r, celldata.connec, ...
            celldata.EdgeData, celldata.verttocell, param, T1flagVec, T1relaxstepcountVec, 0);
    end
    
    Coordinates(:, 2*tstep - 1:2*tstep) = celldata.r;
    Connectivity(:, tstep) = celldata.connec;
    
    param.Tsim = param.Tsim + param.deltat;
    celldata.r0 = celldata.r;
    
    energymat(tstep) = getTissueEnergyClassical(celldata, param);
    timemat(tstep) = param.Tsim;
    if tstep == 1
        relenergychange = (energymat(tstep) - energymat0) / energymat0;
    else
        relenergychange = (energymat(tstep) - energymat(tstep-1)) / energymat(tstep-1);
    end
    
    PrintInformation;
    
    %% Compute and store edge forces at snapshot timesteps
    if ismember(tstep, snapshotTimesteps)
        snapshotCount = snapshotCount + 1;
        edgeForceSnapshots(snapshotCount).tstep = tstep;
        edgeForceSnapshots(snapshotCount).data = getEdgeForces(celldata, param, tstep);
        fprintf('Edge force snapshot saved at timestep %d\n', tstep);
    end
end

simTime = toc(totalsimtime);
fprintf('Total simulation time: %.2f minutes\n', simTime / 60);

%% 4 - Save results
saveDir = 'results/edge_forces';
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end

% --- Plot 1: Final state edge force distribution (4-panel) ---
fig1 = figure('visible', 'off', 'Position', [100 100 1400 1000]);
finalEdgeForces = getEdgeForces(celldata, param, param.Nsteps);
plotEdgeForces(celldata, param, finalEdgeForces, param.Nsteps, 'all');
sgtitle(sprintf('Edge Force Distribution at Final Timestep (%d)', param.Nsteps), 'FontSize', 14);
saveas(fig1, fullfile(saveDir, 'edge_forces_final_4panel.png'));
close(fig1);

% --- Plot 2: Force vectors at final state ---
fig2 = figure('visible', 'off', 'Position', [100 100 800 700]);
plotEdgeForces(celldata, param, finalEdgeForces, param.Nsteps, 'vectors');
title(sprintf('Edge Force Vectors at Timestep %d', param.Nsteps));
saveas(fig2, fullfile(saveDir, 'edge_force_vectors_final.png'));
close(fig2);

% --- Plot 3: Edge tension evolution across snapshots ---
fig3 = figure('visible', 'off', 'Position', [100 100 1400 300*snapshotCount]);
for s = 1:snapshotCount
    subplot(snapshotCount, 1, s);
    % Reconstruct celldata at this snapshot from stored Coordinates/Connectivity
    ts = edgeForceSnapshots(s).tstep;
    tmpCelldata = celldata;
    tmpCelldata.r = Coordinates(:, 2*ts-1:2*ts);
    tmpCelldata.connec = Connectivity(:, ts);
    plotEdgeColorMap_standalone(tmpCelldata, param, edgeForceSnapshots(s).data, 'total_tension');
    title(sprintf('Total Edge Tension at t=%d', ts));
end
sgtitle('Edge Tension Evolution Over Time', 'FontSize', 14);
saveas(fig3, fullfile(saveDir, 'edge_tension_evolution.png'));
close(fig3);

% --- Plot 4: Energy vs time ---
fig4 = figure('visible', 'off');
plot(timemat, energymat, 'LineWidth', 3.0);
title('Energy vs Time');
xlabel('Time');
ylabel('Free Energy');
saveas(fig4, fullfile(saveDir, 'energy_vs_time.png'));
close(fig4);

% --- Save workspace ---
save(fullfile(saveDir, 'edge_force_workspace.mat'));

fprintf('\n=== All edge force results saved to %s ===\n', saveDir);


%% ========================================================================
% Local helper for standalone edge coloring (avoids nested function issues)
% =========================================================================
function plotEdgeColorMap_standalone(celldata, param, edgeForceData, fieldName)

Lx = param.Lx;
Ly = param.Ly;
nEdges = length(edgeForceData);

values = zeros(nEdges, 1);
for i = 1:nEdges
    values(i) = edgeForceData(i).(fieldName);
end

validIdx = [edgeForceData.length] > 1e-10;
hold on;
rectangle('Position', [0 0 Lx Ly], 'EdgeColor', [0.5 0.5 0.5], 'LineStyle', '--');

validValues = values(validIdx);
if isempty(validValues) || (max(validValues) - min(validValues)) < 1e-12
    cmin = -1; cmax = 1;
else
    cmax = max(abs(validValues));
    cmin = -cmax;
end

for i = 1:nEdges
    if ~validIdx(i), continue; end
    r1 = celldata.r(edgeForceData(i).v1, :);
    eVec = edgeForceData(i).edgeVec;
    r2_local = r1 + eVec;
    if abs(eVec(1)) > Lx/2 || abs(eVec(2)) > Ly/2, continue; end
    
    t = (values(i) - cmin) / (cmax - cmin);
    t = max(0, min(1, t));
    if t < 0.5
        s = t / 0.5;
        c = [s, s, 1];
    else
        s = (t - 0.5) / 0.5;
        c = [1, 1-s, 1-s];
    end
    
    lw = 1.0 + 3.0 * abs(values(i)) / max(abs(cmax), 1e-6);
    lw = min(lw, 5.0);
    plot([r1(1) r2_local(1)], [r1(2) r2_local(2)], '-', 'Color', c, 'LineWidth', lw);
    
    if edgeForceData(i).isWoundEdge
        mp = edgeForceData(i).midpoint;
        plot(mp(1), mp(2), 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'k');
    end
end

half = 128;
r1c = linspace(0,1,half)'; g1c = linspace(0,1,half)'; b1c = ones(half,1);
r2c = ones(half,1); g2c = linspace(1,0,half)'; b2c = linspace(1,0,half)';
cmap = [r1c g1c b1c; r2c g2c b2c];
colormap(cmap);
caxis([cmin cmax]);
cb = colorbar;
cb.Label.String = strrep(fieldName, '_', ' ');
axis equal;
xlim([-0.5 Lx+0.5]); ylim([-0.5 Ly+0.5]);
hold off;
end
