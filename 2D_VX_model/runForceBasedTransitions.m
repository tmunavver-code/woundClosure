%% runForceBasedTransitions: Force-Based Transitions — Comprehensive Results
%
% Follows the strategy in README_force_based_transitions.md:
%   Phase 1: Instrument & validate the force field
%   Phase 2: Force-based transition criterion (Bell's law)
%
% Saves ALL results (plots, snapshots, workspace) to:
%   results/wound_contractility/
%
% Key outputs:
%   - Tissue snapshot photos at multiple timesteps (with forces per length)
%   - Force-vs-length scatter at T1 transitions
%   - Trepat spatial force profile (time-averaged)
%   - Edge tension evolution color maps
%   - Energy vs time
%   - Per-edge force-per-length histograms
%   - Wound area evolution
%
% Usage: Run in MATLAB from the 2D_VX_model directory.

addpath(genpath('../common'))
addpath(genpath('./functions2D'))

close all;
clear;
clc;

rng('default');

%% ========================================================================
%  1 — PARAMETERS
% =========================================================================
param.Lx = 15;
param.Ly = 15;

Case = 'fluid';

SingleCellContractionParameters;

% --- Override: Enable T1 transitions for force-based testing ---
param.enableT1transitions = true;
param.enhanceT1atWound    = true;

% --- TEMP: shortened run for quick validation (revert to 3000 for full runs) ---
param.Nsteps = 1000;

% --- Force-based transitions are ON (set in SingleCellContractionParameters) ---
fprintf('=== Force-Based Transitions: %s ===\n', ...
    iif(param.useForceBasedTransitions, 'ENABLED', 'DISABLED'));
fprintf('   k_off0 = %.3f, f_beta = %.3f, f_beta_length_slope = %.3f\n', ...
    param.k_off0, param.f_beta, param.f_beta_length_slope);
fprintf('   Case = %s, p0 = %.2f\n', Case, param.p0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% --- WOUND CONFIGURATION --- %%%
param.woundType = 'circular';
param.woundCircular.center = [param.Lx/2, param.Ly/2];
param.woundCircular.radius = 2.0;
param.woundElliptical.center = [param.Lx/2, param.Ly/2];
param.woundElliptical.x_axis = 8.0;
param.woundElliptical.y_axis = 1.0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ========================================================================
%  2 — OUTPUT DIRECTORY
% =========================================================================
saveDir = 'results/wound_contractility';
if ~exist(saveDir, 'dir')
    mkdir(saveDir);
end
fprintf('Results will be saved to: %s\n', saveDir);

%% ========================================================================
%  3 — GENERATE MESH & CREATE WOUND
% =========================================================================
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

%% ========================================================================
%  4 — RUN SIMULATION
% =========================================================================
totalsimtime = tic;

energymat = zeros(param.Nsteps, 1);
timemat   = zeros(param.Nsteps, 1);
T1flagVec = zeros(celldata.nMasterVertices, 1);
T1relaxstepcountVec = zeros(size(T1flagVec));

Coordinates  = zeros(celldata.nMasterVertices, 2*param.Nsteps);
Connectivity = cell(celldata.nCells, param.Nsteps);

% --- Snapshot schedule: many more snapshots for rich visualization ---
snapshotInterval = max(1, round(param.Nsteps / 20));   % ~20 snapshots
snapshotTimesteps = unique([1, snapshotInterval:snapshotInterval:param.Nsteps, param.Nsteps]);
edgeForceSnapshots = struct();
snapshotCount = 0;

% --- T1 Force logging ---
T1ForceLog = [];

% --- Track wound area over time ---
woundAreaLog = zeros(param.Nsteps, 1);

% --- Track per-step force-per-length statistics ---
meanTangPerLength = zeros(param.Nsteps, 1);
meanNormPerLength = zeros(param.Nsteps, 1);
maxTangPerLength  = zeros(param.Nsteps, 1);
maxNormPerLength  = zeros(param.Nsteps, 1);

fprintf('\n--- Starting simulation: %d timesteps ---\n', param.Nsteps);

for tstep = 1:param.Nsteps
    
    applyStretchX;
    
    [energymat0] = getTissueEnergyClassical(celldata, param);
    
    celldata.f = getVertexForcesClassical(celldata, param, tstep);
    
    % --- Compute edge forces ONCE per step (Phase 1b fix) ---
    currentEdgeForces = getEdgeForces(celldata, param, tstep);
    
    % --- Collect force-per-length stats ---
    validEdges = [currentEdgeForces.length] > 1e-10;
    if any(validEdges)
        tangMags = arrayfun(@(e) norm(e.force_tangent_per_length), currentEdgeForces(validEdges));
        normMags = arrayfun(@(e) norm(e.force_normal_per_length),  currentEdgeForces(validEdges));
        meanTangPerLength(tstep) = mean(tangMags);
        meanNormPerLength(tstep) = mean(normMags);
        maxTangPerLength(tstep)  = max(tangMags);
        maxNormPerLength(tstep)  = max(normMags);
    end
    
    % --- Wound intercalations (T2) ---
    if param.enableWoundIntercalations
        [celldata, T1flagVec, intercalation_happened, param] = ...
            checkWoundIntercalations(celldata, param, T1flagVec, currentEdgeForces, tstep);
        if intercalation_happened
            celldata.A = getCellAreas(celldata, param);
            [celldata.P, celldata.EdgeData] = getCellPerimeters(...
                celldata.nCells, celldata.r, celldata.connec, param, 0);
        end
    end
    
    celldata = updateVertexPositions(celldata, param);
    
    celldata.A = getCellAreas(celldata, param);
    [celldata.P, celldata.EdgeData] = getCellPerimeters(...
        celldata.nCells, celldata.r, celldata.connec, param, 0);
    
    % --- T1 transitions (force-based) ---
    if param.enableT1transitions
        [celldata.r, celldata.connec, celldata.EdgeData, celldata.verttocell, ...
         T1flagVec, T1relaxstepcountVec, param.nT1, T1ForceLog] = ...
            checkT1transitions(celldata.nCells, celldata.r, celldata.connec, ...
            celldata.EdgeData, celldata.verttocell, param, T1flagVec, ...
            T1relaxstepcountVec, 0, currentEdgeForces, tstep, T1ForceLog, celldata);
    end
    
    Coordinates(:, 2*tstep - 1:2*tstep) = celldata.r;
    Connectivity(:, tstep) = celldata.connec;
    
    param.Tsim = param.Tsim + param.deltat;
    celldata.r0 = celldata.r;
    
    energymat(tstep) = getTissueEnergyClassical(celldata, param);
    timemat(tstep)   = param.Tsim;
    if tstep == 1
        relenergychange = (energymat(tstep) - energymat0) / energymat0;
    else
        relenergychange = (energymat(tstep) - energymat(tstep-1)) / energymat(tstep-1);
    end
    
    % --- Wound area: trace boundary as a connected edge chain ---
    if isfield(celldata, 'woundEdges') && ~isempty(celldata.woundEdges)
        woundEdgesMat = celldata.woundEdges;
        nWE = size(woundEdgesMat, 1);
        if nWE >= 3
            % Build adjacency from wound edges and trace the loop
            % woundEdges rows are [v1, v2] free-edge pairs
            try
                visited = false(nWE, 1);
                orderedVerts = woundEdgesMat(1, 1);
                nextV = woundEdgesMat(1, 2);
                visited(1) = true;
                maxIter = nWE + 1;
                iter = 0;
                while nextV ~= orderedVerts(1) && iter < maxIter
                    iter = iter + 1;
                    orderedVerts(end+1) = nextV; %#ok<AGROW>
                    % Find the next unvisited edge connected to nextV
                    found = false;
                    for ei = 1:nWE
                        if visited(ei), continue; end
                        if woundEdgesMat(ei,1) == nextV
                            visited(ei) = true;
                            nextV = woundEdgesMat(ei,2);
                            found = true; break;
                        elseif woundEdgesMat(ei,2) == nextV
                            visited(ei) = true;
                            nextV = woundEdgesMat(ei,1);
                            found = true; break;
                        end
                    end
                    if ~found, break; end
                end
                if length(orderedVerts) >= 3
                    px = celldata.r(orderedVerts, 1);
                    py = celldata.r(orderedVerts, 2);
                    % Periodic wrap correction for each coordinate
                    px = px - param.Lx * round((px - px(1)) / param.Lx);
                    py = py - param.Ly * round((py - py(1)) / param.Ly);
                    woundAreaLog(tstep) = abs(polyarea(px, py));
                end
            catch
                % Fall back to convex hull area if chain tracing fails
                woundVerts = unique(woundEdgesMat(:));
                woundCoords = celldata.r(woundVerts, :);
                if size(woundCoords, 1) >= 3
                    try
                        k = convhull(woundCoords(:,1), woundCoords(:,2));
                        woundAreaLog(tstep) = polyarea(woundCoords(k,1), woundCoords(k,2));
                    catch
                        woundAreaLog(tstep) = woundAreaLog(max(1,tstep-1));
                    end
                end
            end
        end
    end
    
    PrintInformation;
    
    %% --- Save snapshots ---
    if ismember(tstep, snapshotTimesteps)
        snapshotCount = snapshotCount + 1;
        edgeForceSnapshots(snapshotCount).tstep = tstep;
        edgeForceSnapshots(snapshotCount).data  = currentEdgeForces;
        
        % ============================================================
        %  SNAPSHOT PHOTO: tissue + forces per length at this timestep
        % ============================================================
        figSnap = figure('visible', 'off', 'Position', [50 50 1200 500]);
        
        % Panel 1: Tissue with force-per-length vectors
        subplot(1, 2, 1);
        plotEdgeForces(celldata, param, currentEdgeForces, tstep, 'vectors_per_length');
        title(sprintf('Force / Length  (t = %d)', tstep), 'FontSize', 12);
        
        % Panel 2: Edge tension colormap
        subplot(1, 2, 2);
        plotEdgeColorMap_standalone(celldata, param, currentEdgeForces, 'total_tension');
        title(sprintf('Total Tension  (t = %d)', tstep), 'FontSize', 12);
        
        sgtitle(sprintf('Timestep %d / %d  |  Force-Based Transitions', tstep, param.Nsteps), ...
            'FontSize', 14, 'FontWeight', 'bold');
        
        saveas(figSnap, fullfile(saveDir, sprintf('snapshot_t%04d.png', tstep)));
        close(figSnap);
        
        fprintf('  >> Snapshot saved: snapshot_t%04d.png\n', tstep);
    end
end

simTime = toc(totalsimtime);
fprintf('\n=== Simulation finished in %.2f minutes ===\n', simTime / 60);
fprintf('Total T1 transitions: %d\n', param.nT1);

%% ========================================================================
%  5 — POST-PROCESSING & RESULT FIGURES
% =========================================================================

finalEdgeForces = getEdgeForces(celldata, param, param.Nsteps);

% --- Fig 1: Final state 4-panel force distribution ---
fig1 = figure('visible', 'off', 'Position', [100 100 1400 1000]);
plotEdgeForces(celldata, param, finalEdgeForces, param.Nsteps, 'all');
sgtitle(sprintf('Edge Force Distribution — Final (t=%d) | p0=%.2f', ...
    param.Nsteps, param.p0), 'FontSize', 14);
saveas(fig1, fullfile(saveDir, '01_edge_forces_final_4panel.png'));
close(fig1);

% --- Fig 2: Force vectors per unit length at final state ---
fig2 = figure('visible', 'off', 'Position', [100 100 900 800]);
plotEdgeForces(celldata, param, finalEdgeForces, param.Nsteps, 'vectors_per_length');
title(sprintf('Edge Force Vectors / Length — Final (t=%d)', param.Nsteps));
saveas(fig2, fullfile(saveDir, '02_force_per_length_final.png'));
close(fig2);

% --- Fig 3: Force vectors (total) at final state ---
fig3 = figure('visible', 'off', 'Position', [100 100 900 800]);
plotEdgeForces(celldata, param, finalEdgeForces, param.Nsteps, 'vectors');
title(sprintf('Edge Force Vectors — Final (t=%d)', param.Nsteps));
saveas(fig3, fullfile(saveDir, '03_force_vectors_final.png'));
close(fig3);

% --- Fig 4: Edge tension evolution across snapshots ---
nSnaps = min(snapshotCount, 10); % cap at 10 rows for readability
fig4 = figure('visible', 'off', 'Position', [100 100 1400  250*nSnaps]);
for s = 1:nSnaps
    subplot(nSnaps, 1, s);
    ts = edgeForceSnapshots(s).tstep;
    tmpCelldata = celldata;
    tmpCelldata.r = Coordinates(:, 2*ts-1:2*ts);
    tmpCelldata.connec = Connectivity(:, ts);
    plotEdgeColorMap_standalone(tmpCelldata, param, edgeForceSnapshots(s).data, 'total_tension');
    title(sprintf('Total Edge Tension at t = %d', ts));
end
sgtitle('Edge Tension Evolution Over Time', 'FontSize', 14);
saveas(fig4, fullfile(saveDir, '04_tension_evolution.png'));
close(fig4);

% --- Fig 5: Energy vs time ---
fig5 = figure('visible', 'off', 'Position', [100 100 700 400]);
plot(timemat, energymat, 'LineWidth', 2.0, 'Color', [0.2 0.4 0.8]);
title('Free Energy vs Time'); xlabel('Time'); ylabel('Free Energy');
grid on;
saveas(fig5, fullfile(saveDir, '05_energy_vs_time.png'));
close(fig5);

% --- Fig 6: T1 Force-vs-Length scatter (Phase 1 diagnostic) ---
if ~isempty(T1ForceLog)
    fig6 = figure('visible', 'off', 'Position', [100 100 1100 450]);
    subplot(1,2,1);
    scatter(T1ForceLog(:,4), T1ForceLog(:,5), 50, T1ForceLog(:,1), 'filled');
    xlabel('Edge Length at Flip'); ylabel('Tangential Tension / Length');
    title('T1: Tangential Force / Length'); grid on; colorbar;
    cb = colorbar; cb.Label.String = 'Timestep';
    
    subplot(1,2,2);
    scatter(T1ForceLog(:,4), T1ForceLog(:,6), 50, T1ForceLog(:,1), 'filled');
    xlabel('Edge Length at Flip'); ylabel('Normal Pressure');
    title('T1: Normal Pressure'); grid on;
    cb = colorbar; cb.Label.String = 'Timestep';
    
    sgtitle('Forces at the Moment of T1 Transition (Phase 1 Diagnostic)');
    saveas(fig6, fullfile(saveDir, '06_t1_forces_vs_length.png'));
    close(fig6);
else
    fprintf('No T1 transitions occurred — skipping force-vs-length scatter.\n');
end

% --- Fig 7: Trepat Validation — single snapshot ---
fig7 = validateTrepatForces(finalEdgeForces, param, param.Nsteps);
saveas(fig7, fullfile(saveDir, '07_trepat_single_snapshot.png'));
close(fig7);

% --- Fig 8: Trepat Validation — TIME-AVERAGED (Phase 1 requirement) ---
fig8 = validateTrepatForces(edgeForceSnapshots, param, param.Nsteps, 'timeAverage', true);
saveas(fig8, fullfile(saveDir, '08_trepat_time_averaged.png'));
close(fig8);

% --- Fig 9: Force-per-length statistics over time ---
fig9 = figure('visible', 'off', 'Position', [100 100 1000 500]);
subplot(2,1,1);
plot(timemat, meanTangPerLength, 'g-', 'LineWidth', 1.5); hold on;
plot(timemat, maxTangPerLength,  'g--', 'LineWidth', 1.0);
ylabel('Tangential F/L'); title('Force per Length Statistics over Time');
legend('Mean |F_{tan}/L|', 'Max |F_{tan}/L|', 'Location', 'best');
grid on; hold off;

subplot(2,1,2);
plot(timemat, meanNormPerLength, 'b-', 'LineWidth', 1.5); hold on;
plot(timemat, maxNormPerLength,  'b--', 'LineWidth', 1.0);
ylabel('Normal F/L'); xlabel('Time');
legend('Mean |F_{norm}/L|', 'Max |F_{norm}/L|', 'Location', 'best');
grid on; hold off;

saveas(fig9, fullfile(saveDir, '09_force_per_length_timeseries.png'));
close(fig9);

% --- Fig 10: Force-per-length histogram at final timestep ---
fig10 = figure('visible', 'off', 'Position', [100 100 1000 400]);
validFinal = [finalEdgeForces.length] > 1e-10;
tangMagsFinal = arrayfun(@(e) norm(e.force_tangent_per_length), finalEdgeForces(validFinal));
normMagsFinal = arrayfun(@(e) norm(e.force_normal_per_length),  finalEdgeForces(validFinal));

subplot(1,2,1);
histogram(tangMagsFinal, 30, 'FaceColor', [0.2 0.8 0.2], 'EdgeColor', 'k');
xlabel('|F_{tan}| / L'); ylabel('Count');
title('Tangential Force / Length Distribution');
grid on;

subplot(1,2,2);
histogram(normMagsFinal, 30, 'FaceColor', [0.2 0.4 0.9], 'EdgeColor', 'k');
xlabel('|F_{norm}| / L'); ylabel('Count');
title('Normal Force / Length Distribution');
grid on;

sgtitle(sprintf('Force per Length Histograms — Final (t=%d)', param.Nsteps));
saveas(fig10, fullfile(saveDir, '10_force_per_length_histogram.png'));
close(fig10);

% --- Fig 11: Wound area evolution ---
fig11 = figure('visible', 'off', 'Position', [100 100 700 400]);
plot(timemat, woundAreaLog, 'r-', 'LineWidth', 2.0);
title('Wound Area vs Time'); xlabel('Time'); ylabel('Wound Area');
grid on;
saveas(fig11, fullfile(saveDir, '11_wound_area_vs_time.png'));
close(fig11);

% --- Fig 12: Wound-edge vs bulk force-per-length comparison ---
fig12 = figure('visible', 'off', 'Position', [100 100 900 400]);
woundIdx  = [finalEdgeForces.isWoundEdge] & validFinal;
bulkIdx   = ~[finalEdgeForces.isWoundEdge] & validFinal;

tangWound = arrayfun(@(e) norm(e.force_tangent_per_length), finalEdgeForces(woundIdx));
tangBulk  = arrayfun(@(e) norm(e.force_tangent_per_length), finalEdgeForces(bulkIdx));
normWound = arrayfun(@(e) norm(e.force_normal_per_length),  finalEdgeForces(woundIdx));
normBulk  = arrayfun(@(e) norm(e.force_normal_per_length),  finalEdgeForces(bulkIdx));

subplot(1,2,1);
barData = [mean(tangWound), mean(tangBulk)];
barErr  = [std(tangWound),  std(tangBulk)];
b = bar(barData, 'FaceColor', 'flat'); hold on;
b.CData(1,:) = [0.9 0.3 0.3]; b.CData(2,:) = [0.3 0.6 0.9];
errorbar(1:2, barData, barErr, 'k.', 'LineWidth', 1.5);
set(gca, 'XTickLabel', {'Wound Edge', 'Bulk'});
ylabel('|F_{tan}| / L'); title('Tangential F/L: Wound vs Bulk');
grid on; hold off;

subplot(1,2,2);
barData2 = [mean(normWound), mean(normBulk)];
barErr2  = [std(normWound),  std(normBulk)];
b2 = bar(barData2, 'FaceColor', 'flat'); hold on;
b2.CData(1,:) = [0.9 0.3 0.3]; b2.CData(2,:) = [0.3 0.6 0.9];
errorbar(1:2, barData2, barErr2, 'k.', 'LineWidth', 1.5);
set(gca, 'XTickLabel', {'Wound Edge', 'Bulk'});
ylabel('|F_{norm}| / L'); title('Normal F/L: Wound vs Bulk');
grid on; hold off;

sgtitle('Wound-Edge vs Bulk Force per Length Comparison');
saveas(fig12, fullfile(saveDir, '12_wound_vs_bulk_force_per_length.png'));
close(fig12);

%% ========================================================================
%  6 — SAVE WORKSPACE & SUMMARY
% =========================================================================
save(fullfile(saveDir, 'workspace.mat'));

% Write a summary text file
fid = fopen(fullfile(saveDir, 'RUN_SUMMARY.txt'), 'w');
fprintf(fid, 'Force-Based Transitions Run — Comprehensive Results\n');
fprintf(fid, '==========================================================\n');
fprintf(fid, 'Date: %s\n', datestr(now));
fprintf(fid, 'Case: %s\n', Case);
fprintf(fid, 'p0 = %.2f\n', param.p0);
fprintf(fid, 'Lx x Ly = %.0f x %.0f\n', param.Lx, param.Ly);
fprintf(fid, 'Nsteps = %d, dt = %.4f\n', param.Nsteps, param.deltat);
fprintf(fid, 'useForceBasedTransitions = %d\n', param.useForceBasedTransitions);
fprintf(fid, 'k_off0 = %.3f, f_beta = %.3f, f_beta_length_slope = %.3f\n', ...
    param.k_off0, param.f_beta, param.f_beta_length_slope);
fprintf(fid, 'enableT1transitions = %d\n', param.enableT1transitions);
fprintf(fid, 'enableWoundIntercalations = %d\n', param.enableWoundIntercalations);
fprintf(fid, 'lambda_purse_string = %.2f\n', param.lambda_purse_string);
fprintf(fid, 'T1_TOL_bulk = %.3f, T1_TOL_wound = %.3f\n', param.T1_TOL_bulk, param.T1_TOL_wound);
if isfield(param, 'ka_wound_factor')
    fprintf(fid, 'ka_wound_factor = %.2f\n', param.ka_wound_factor);
end
if isfield(param, 'contractility_wound')
    fprintf(fid, 'contractility_wound = %.2f\n', param.contractility_wound);
end
if isfield(param, 'enableWoundCrawling')
    fprintf(fid, 'enableWoundCrawling = %d, crawl_force0 = %.2f, tau_decay_crawl = %.1f\n', ...
        param.enableWoundCrawling, param.crawl_force0, param.tau_decay_crawl);
end
fprintf(fid, '\nRESULTS:\n');
fprintf(fid, 'Total T1 transitions: %d\n', param.nT1);
fprintf(fid, 'Total simulation time: %.2f minutes\n', simTime / 60);
fprintf(fid, 'Final energy: %.4f\n', energymat(end));
fprintf(fid, 'Final wound area: %.4f\n', woundAreaLog(end));
fprintf(fid, '\nSnapshots saved at timesteps: %s\n', mat2str(snapshotTimesteps));
fprintf(fid, '\nFILES IN THIS DIRECTORY:\n');
fprintf(fid, '  01_edge_forces_final_4panel.png\n');
fprintf(fid, '  02_force_per_length_final.png\n');
fprintf(fid, '  03_force_vectors_final.png\n');
fprintf(fid, '  04_tension_evolution.png\n');
fprintf(fid, '  05_energy_vs_time.png\n');
fprintf(fid, '  06_t1_forces_vs_length.png\n');
fprintf(fid, '  07_trepat_single_snapshot.png\n');
fprintf(fid, '  08_trepat_time_averaged.png\n');
fprintf(fid, '  09_force_per_length_timeseries.png\n');
fprintf(fid, '  10_force_per_length_histogram.png\n');
fprintf(fid, '  11_wound_area_vs_time.png\n');
fprintf(fid, '  12_wound_vs_bulk_force_per_length.png\n');
fprintf(fid, '  snapshot_tXXXX.png (tissue + force/length at each snapshot)\n');
fprintf(fid, '  workspace.mat\n');
fclose(fid);

fprintf('\n=== All results saved to %s ===\n', saveDir);
fprintf('Files saved:\n');
dir(fullfile(saveDir, '*.png'));


%% ========================================================================
% LOCAL HELPERS
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

function result = iif(condition, trueVal, falseVal)
if condition
    result = trueVal;
else
    result = falseVal;
end
end
