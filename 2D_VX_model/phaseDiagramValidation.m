%% phaseDiagramValidation: reproduce the Farhadifar et al. (2007) ground-state
%  phase diagram with THIS repository's own mechanics functions.
%
%  Farhadifar energy per cell:  K/2 (A - A0)^2 + Gamma/2 P^2 + Lambda * sum(edges)
%  This model's energy per cell: ka (A - 1)^2 + (1/rstiff) (P - p0)^2
%
%  Expanding and matching term by term (A0 = 1, each edge shared by two cells):
%     K = 2 ka,   Gamma = 2/rstiff,   Lambda = -4 p0 / rstiff
%  so the dimensionless axes of the figure are
%     Gbar = Gamma / (K A0)           = 1 / (ka * rstiff)
%     Lbar = Lambda / (K A0^(3/2))    = -2 * p0 / (ka * rstiff) = -2 * p0 * Gbar
%  and the soft/hexagonal boundary is the ray p0 = sqrt(8*sqrt(3)) = 3.7224,
%  i.e. Lbar = -2 * 3.7224 * Gbar.
%
%  For each (Lbar, Gbar) point: start from a periodic hexagonal lattice of
%  unit-area cells with a small random perturbation, relax it with
%  getVertexForcesClassical + updateVertexPositions (no wound, no T1s), and
%  classify the relaxed state by its mean shape index q = <P / sqrt(A)>:
%     hexagonal : q stays at 3.7224 (regular hexagons)
%     soft      : q rises above it (irregular, degenerate network)
%     collapse  : some cell area falls below 0.3
%
%  Part 0 is a unit test: getVertexForcesClassical must equal the numerical
%  gradient of getTissueEnergyClassical.
%
%  The box is fixed (as in every simulation in this repository), so the
%  hatched "collapse" region of the published figure, which assumes a
%  pressure-free box, is not tested here. In the fixed box the hexagonal
%  network stays stable well past that boundary; Part 3 checks this up to
%  Lbar = 2.4 (the NumPy port of these functions stays hexagonal to 2.7).
%
%  Outputs go to results/phase_diagram_validation_<timestamp>/.

clear; close all; clc;

% Run from this folder WITHOUT adding functions2D/ to the path: several files
% there are stale copies that would shadow the current versions.
thisDir = fileparts(mfilename('fullpath'));
cd(thisDir);

runMode = 'quick';   % 'quick' (~50 points, well under an hour) or 'full'
dt      = 0.05;      % explicit-Euler step. The model default 0.01 gives the
                     % same phases (identical energy landscape), 5x slower.
nx = 6; ny = 6;      % 36-cell periodic lattice

c = sqrt(8*sqrt(3));
saveDir = fullfile('results', ['phase_diagram_validation_' datestr(now,'yyyy-mm-dd_HHMMSS')]);
mkdir(saveDir);
fprintf('Saving to %s\n', saveDir);

%% Part 0 - force = -grad(energy) unit test on the repository's functions
[r, connec, Lx, Ly] = buildHexLattice(nx, ny);
rng(0);
r = r + 0.05*randn(size(r));
r(:,1) = mod(r(:,1), Lx); r(:,2) = mod(r(:,2), Ly);
param = makeParam(Lx, Ly, -0.75, 0.10, dt);
cd0 = makeCelldata(r, connec, param);
f = getVertexForcesClassical(cd0, param, 1);
h = 1e-6; g = zeros(size(r));
for i = 1:size(r,1)
    for k = 1:2
        rp = r; rp(i,k) = rp(i,k) + h;
        rm = r; rm(i,k) = rm(i,k) - h;
        g(i,k) = (getTissueEnergyClassical(makeCelldata(rp,connec,param),param) - ...
                  getTissueEnergyClassical(makeCelldata(rm,connec,param),param)) / (2*h);
    end
end
forceErr = max(abs(f(:) + g(:))) / max(abs(g(:)));
fprintf('Part 0: max relative |f + dE/dr| = %.3e  (should be ~1e-9)\n', forceErr);

%% Point list
pts = struct('set', {}, 'Lbar', {}, 'Gbar', {}, 'nmax', {});
switch runMode
    case 'quick'
        lineG  = [0.05 0.11 0.17];
        lineP0 = [3.62 3.69 3.71 3.72 3.725 3.735 3.76 3.82];
        gridG  = [0.03 0.09 0.15];
        gridL  = -1.4:0.4:0.6;
    otherwise
        lineG  = [0.02 0.05 0.08 0.11 0.14 0.17];
        lineP0 = [3.55 3.62 3.67 3.70 3.71 3.715 3.72 3.725 3.73 3.74 3.76 3.80 3.88];
        gridG  = 0.01:0.02:0.17;
        gridL  = -1.4:0.1:0.6;
end
for G = lineG
    for p0 = lineP0
        pts(end+1) = struct('set','line','Lbar',-2*p0*G,'Gbar',G,'nmax',round(800/dt)); %#ok<SAGROW>
    end
end
labelled = {'I',0.12,0.04; 'II',0.0,0.10; 'III',-0.85,0.10; 'IV',-0.32,0.04; 'V',0.0,0.04};
for k = 1:size(labelled,1)
    pts(end+1) = struct('set',['point_' labelled{k,1}],'Lbar',labelled{k,2}, ...
        'Gbar',labelled{k,3},'nmax',round(800/dt)); %#ok<SAGROW>
end
for G = gridG
    for L = gridL
        pts(end+1) = struct('set','grid','Lbar',L,'Gbar',G,'nmax',round(300/dt)); %#ok<SAGROW>
    end
end
for G = [0.05 0.12]
    for L = [1.9 2.05 2.25 2.4]
        pts(end+1) = struct('set','fixedbox_collapse','Lbar',L,'Gbar',G,'nmax',round(400/dt)); %#ok<SAGROW>
    end
end

%% Parts 1-3 - relax every point
nP = numel(pts);
res = table('Size',[nP 12], 'VariableTypes', ...
    {'string','double','double','double','double','string','double','double','double','double','double','string'}, ...
    'VariableNames', {'set','Lbar','Gbar','p0','steps','status','q_mean','q_minus_c', ...
    'edge_cv','energy_per_cell','A_min','phase'});
tAll = tic;
for n = 1:nP
    [row, cdFinal] = relaxPoint(pts(n).Lbar, pts(n).Gbar, nx, ny, pts(n).nmax, dt, c);
    res(n,:) = {pts(n).set, pts(n).Lbar, pts(n).Gbar, row.p0, row.steps, row.status, ...
        row.q, row.q - c, row.cv, row.epc, row.Amin, row.phase};
    fprintf('[%3d/%3d] %-18s Lbar=%+.3f Gbar=%.3f p0=%+.3f -> %-9s q-c=%+.2e (%s, %d steps, %.0f s total)\n', ...
        n, nP, pts(n).set, pts(n).Lbar, pts(n).Gbar, row.p0, row.phase, row.q - c, ...
        row.status, row.steps, toc(tAll));
    if strcmp(pts(n).set, 'point_III'), cdSoft = cdFinal; paramSoft = makeParam(Lx,Ly,pts(n).Lbar,pts(n).Gbar,dt); end
    if strcmp(pts(n).set, 'point_II'),  cdHex  = cdFinal; paramHex  = makeParam(Lx,Ly,pts(n).Lbar,pts(n).Gbar,dt); end
end
writetable(res, fullfile(saveDir, 'phase_diagram_points.csv'));

%% Boundary estimate along the line: midpoint between last hexagonal and first soft p0
fprintf('\nSoft/hexagonal boundary (theory p0* = %.4f):\n', c);
for G = lineG
    sel = res.set == "line" & abs(res.Gbar - G) < 1e-12;
    p = res.p0(sel); ph = res.phase(sel);
    hexMax = max(p(ph == "hexagonal")); softMin = min(p(ph == "soft"));
    fprintf('  Gbar = %.2f : hexagonal up to p0 = %.3f, soft from p0 = %.3f -> p0* in [%.3f, %.3f]\n', ...
        G, hexMax, softMin, hexMax, softMin);
end

%% Figure
fig = figure('Color','white','Position',[80 80 1100 760]);
ax = axes(fig); hold(ax,'on'); box(ax,'on');
Gm = 0.18;
% soft region: left of the p0 = 3.7224 ray
patch(ax, [-1.4 0 -2*c*Gm -1.4], [0 0 Gm Gm], [0.88 0.94 0.98], 'EdgeColor','none');
% hexagonal region, right boundary from the pressure-free hexagonal network
Gc = linspace(0, 2/c^2, 200);
LR = 8/c * ((1 - 0.5*Gc*c^2)/3).^1.5;
patch(ax, [0 -2*c*Gm 0 0 fliplr(LR)], [0 Gm Gm 2/c^2 fliplr(Gc)], [0.80 0.80 0.80], 'EdgeColor','none');
plot(ax, [0 -2*c*Gm], [0 Gm], 'k-', 'LineWidth', 1.5);
plot(ax, [LR 0 0], [Gc 2/c^2 Gm], 'k--', 'LineWidth', 1.2);
isHex  = res.phase == "hexagonal";
isSoft = res.phase == "soft";
onChart = res.Lbar <= 0.6;
scatter(ax, res.Lbar(isHex & onChart),  res.Gbar(isHex & onChart),  28, [0.25 0.25 0.25], 's', 'filled');
scatter(ax, res.Lbar(isSoft & onChart), res.Gbar(isSoft & onChart), 28, [0.10 0.45 0.80], 'o', 'filled');
for k = 1:size(labelled,1)
    plot(ax, labelled{k,2}, labelled{k,3}, 'o', 'MarkerSize', 11, 'MarkerFaceColor', [0.55 0.80 0.30], 'MarkerEdgeColor','k');
    text(ax, labelled{k,2} - 0.03, labelled{k,3} + 0.006, labelled{k,1}, 'FontSize', 12, 'HorizontalAlignment','right');
end
xlim(ax, [-1.4 0.6]); ylim(ax, [0 Gm]);
xlabel(ax, 'normalized line tension  \Lambda / K A_0^{3/2}  =  -2 p_0 / (k_a r_{stiff})');
ylabel(ax, 'normalized contractility  \Gamma / K A_0  =  1 / (k_a r_{stiff})');
title(ax, sprintf('Ground-state phase diagram from this vertex model (force check %.1e)', forceErr));
legend(ax, {'soft (theory)','hexagonal (theory)','p_0 = 3.7224 line', ...
    'collapse boundary (pressure-free box, not tested)','simulated: hexagonal','simulated: soft'}, ...
    'Location','northeast');
exportgraphics(fig, fullfile(saveDir, 'phase_diagram.png'), 'Resolution', 200);
savefig(fig, fullfile(saveDir, 'phase_diagram.fig'));

% relaxed example networks
if exist('cdSoft','var') && exist('cdHex','var')
    fig2 = figure('Color','white','Position',[80 80 1000 480]);
    subplot(1,2,1); plot3Dtissue(cdHex.nCells, cdHex.r, cdHex.connec, paramHex); axis equal off;
    title(sprintf('Point II (\\Lambda=0, \\Gamma=0.10): hexagonal, q = %.4f', mean(cdHex.P./sqrt(cdHex.A))));
    subplot(1,2,2); plot3Dtissue(cdSoft.nCells, cdSoft.r, cdSoft.connec, paramSoft); axis equal off;
    title(sprintf('Point III (\\Lambda=-0.85, \\Gamma=0.10): soft, q = %.4f', mean(cdSoft.P./sqrt(cdSoft.A))));
    exportgraphics(fig2, fullfile(saveDir, 'relaxed_networks.png'), 'Resolution', 200);
end
fprintf('\nDone in %.0f s. Results in %s\n', toc(tAll), saveDir);


%% ---------------------------------------------------------------- helpers
function [r, connec, Lx, Ly] = buildHexLattice(nx, ny)
% Periodic pointy-top hexagonal tiling; every cell has area exactly 1 and
% vertices in anticlockwise order.
a = sqrt(2/(3*sqrt(3)));  w = sqrt(3)*a;  h = 1.5*a;
Lx = nx*w;  Ly = ny*h;
keys = containers.Map('KeyType','char','ValueType','double');
r = zeros(0,2);  connec = cell(nx*ny,1);
ang = deg2rad(30 + 60*(0:5));
cellID = 0;
for j = 0:ny-1
    for i = 0:nx-1
        cellID = cellID + 1;
        cx = (i + 0.5*mod(j,2) + 0.5)*w;  cy = (j + 0.5)*h;
        ids = zeros(1,6);
        for k = 1:6
            vx = mod(cx + a*cos(ang(k)), Lx);  vy = mod(cy + a*sin(ang(k)), Ly);
            key = sprintf('%d_%d', mod(round(vx/Lx*1e6),1e6), mod(round(vy/Ly*1e6),1e6));
            if isKey(keys, key)
                ids(k) = keys(key);
            else
                r(end+1,:) = [vx vy]; %#ok<AGROW>
                keys(key) = size(r,1);
                ids(k) = size(r,1);
            end
        end
        connec{cellID} = ids;
    end
end
end

function param = makeParam(Lx, Ly, Lbar, Gbar, dt)
param.Lx = Lx;  param.Ly = Ly;
param.ka = 1;
param.rstiff = 1/(param.ka*Gbar);
param.p0 = -Lbar/(2*Gbar);
param.deltat = dt;  param.eta = 1;
param.isBoundaryFixed = 0;
param.SelfPropellingCellIDs = 0;
param.cellIDtoContract = [];          % no wound-margin cells
param.cellIDstoTrack = 0;             % read by plot3Dtissue
param.enableWoundCrawling = false;
param.enableWoundResistance = false;
end

function celldata = makeCelldata(r, connec, param)
celldata.nCells = numel(connec);
celldata.nMasterVertices = size(r,1);
celldata.connec = connec;
celldata.r = r;  celldata.r0 = r;  celldata.f = zeros(size(r));
celldata.A = getCellAreas(celldata, param);
[celldata.P, celldata.EdgeData] = getCellPerimeters(celldata.nCells, r, connec, param, 0);
end

function [row, celldata] = relaxPoint(Lbar, Gbar, nx, ny, nmax, dt, c)
[r, connec, Lx, Ly] = buildHexLattice(nx, ny);
rng(1);
r = r + 0.02*randn(size(r));
r(:,1) = mod(r(:,1), Lx);  r(:,2) = mod(r(:,2), Ly);
param = makeParam(Lx, Ly, Lbar, Gbar, dt);
celldata = makeCelldata(r, connec, param);
row.status = 'maxsteps';
for s = 1:nmax
    celldata.f = getVertexForcesClassical(celldata, param, s);
    if any(~isfinite(celldata.f(:))), row.status = 'nonfinite'; break; end
    if max(abs(celldata.f(:))) < 2e-6, row.status = 'converged'; break; end
    celldata = updateVertexPositions(celldata, param);
    celldata.r0 = celldata.r;
    celldata.A = getCellAreas(celldata, param);
    [celldata.P, celldata.EdgeData] = getCellPerimeters(celldata.nCells, celldata.r, celldata.connec, param, 0);
end
row.steps = s;
row.p0 = param.p0;
row.q = mean(celldata.P ./ sqrt(max(celldata.A, 1e-12)));
e = [celldata.EdgeData{:}];
row.cv = std(e)/mean(e);
row.epc = getTissueEnergyClassical(celldata, param) / celldata.nCells;
row.Amin = min(celldata.A);
if row.Amin < 0.3 || strcmp(row.status, 'nonfinite')
    row.phase = 'collapse';
elseif row.q - c > 1e-3
    row.phase = 'soft';
else
    row.phase = 'hexagonal';
end
end
