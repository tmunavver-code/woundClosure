function plotEdgeForces(celldata, param, edgeForceData, tstep, plotType)
%% plotEdgeForces: Visualize force distribution on cell edges
%
% INPUTS:
%   celldata      - cell data struct
%   param         - parameter struct
%   edgeForceData - output from getEdgeForces()
%   tstep         - current timestep (for title)
%   plotType      - string: 'tension', 'pressure', 'total', 'vectors',
%                   'vectors_per_length', or 'all'
%
% PLOT TYPES:
%   'tension'            - Color edges by their line tension (perimeter + purse-string)
%   'pressure'           - Color edges by the area pressure acting across them
%   'total'              - Color edges by total tension magnitude
%   'vectors'            - Draw force arrows (tangential + normal) at edge midpoints
%   'vectors_per_length' - Draw force arrows per unit edge length at edge midpoints
%   'all'                - Generate all four plots in subplots

if nargin < 5
    plotType = 'all';
end

Lx = param.Lx;
Ly = param.Ly;

switch plotType
    case 'all'
        figure('Position', [100 100 1200 900]);
        
        subplot(2,2,1);
        plotEdgeColorMap(celldata, param, edgeForceData, 'tension_peri', tstep);
        title(sprintf('Perimeter Tension (t=%d)', tstep));
        
        subplot(2,2,2);
        plotEdgeColorMap(celldata, param, edgeForceData, 'tension_purse', tstep);
        title(sprintf('Purse-String Tension (t=%d)', tstep));
        
        subplot(2,2,3);
        plotEdgeColorMap(celldata, param, edgeForceData, 'total_tension', tstep);
        title(sprintf('Total Edge Tension (t=%d)', tstep));
        
        subplot(2,2,4);
        plotEdgeColorMap(celldata, param, edgeForceData, 'pressure_area', tstep);
        title(sprintf('Area Pressure (t=%d)', tstep));
        
    case 'tension'
        plotEdgeColorMap(celldata, param, edgeForceData, 'total_tension', tstep);
        title(sprintf('Total Edge Tension (t=%d)', tstep));
        
    case 'pressure'
        plotEdgeColorMap(celldata, param, edgeForceData, 'pressure_area', tstep);
        title(sprintf('Area Pressure (t=%d)', tstep));
        
    case 'vectors'
        plotEdgeVectors(celldata, param, edgeForceData, tstep);
        title(sprintf('Edge Force Vectors (t=%d)', tstep));
        
    case 'vectors_per_length'
        plotEdgeVectorsPerLength(celldata, param, edgeForceData, tstep);
        title(sprintf('Edge Force Vectors per Unit Length (t=%d)', tstep));
        
    otherwise
        plotEdgeColorMap(celldata, param, edgeForceData, 'total_tension', tstep);
        title(sprintf('Total Edge Tension (t=%d)', tstep));
end

end


%% ========================================================================
% HELPER: Color-coded edge map
% =========================================================================
function plotEdgeColorMap(celldata, param, edgeForceData, fieldName, tstep)

Lx = param.Lx;
Ly = param.Ly;
nEdges = length(edgeForceData);

% Extract the scalar field values
values = zeros(nEdges, 1);
for i = 1:nEdges
    values(i) = edgeForceData(i).(fieldName);
end

% Filter out edges with zero length (degenerate)
validIdx = [edgeForceData.length] > 1e-10;

hold on;

% Draw tissue outline
rectangle('Position', [0 0 Lx Ly], 'EdgeColor', [0.5 0.5 0.5], 'LineStyle', '--');

% Determine color range
validValues = values(validIdx);
if isempty(validValues) || (max(validValues) - min(validValues)) < 1e-12
    cmin = -1;
    cmax = 1;
else
    cmax = max(abs(validValues));
    cmin = -cmax; % Symmetric colormap centered on 0
end

% Draw edges colored by value
for i = 1:nEdges
    if ~validIdx(i)
        continue;
    end
    
    v1 = edgeForceData(i).v1;
    v2 = edgeForceData(i).v2;
    r1 = celldata.r(v1, :);
    eVec = edgeForceData(i).edgeVec;
    r2_local = r1 + eVec;
    
    % Skip edges that cross the periodic boundary (would draw long lines)
    if abs(eVec(1)) > Lx/2 || abs(eVec(2)) > Ly/2
        continue;
    end
    
    % Map value to color
    t = (values(i) - cmin) / (cmax - cmin);
    t = max(0, min(1, t));
    
    % Blue (negative/compressive) -> White (zero) -> Red (positive/contractile)
    if t < 0.5
        s = t / 0.5;
        c = [s, s, 1]; % Blue to white
    else
        s = (t - 0.5) / 0.5;
        c = [1, 1 - s, 1 - s]; % White to red
    end
    
    lineWidth = 1.0 + 3.0 * abs(values(i)) / max(abs(cmax), 1e-6);
    lineWidth = min(lineWidth, 5.0);
    
    plot([r1(1) r2_local(1)], [r1(2) r2_local(2)], '-', ...
        'Color', c, 'LineWidth', lineWidth);
    
    % Mark wound edges with a dot
    if edgeForceData(i).isWoundEdge
        mp = edgeForceData(i).midpoint;
        plot(mp(1), mp(2), 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'k');
    end
end

% Add colorbar
colormap(blueWhiteRed(256));
caxis([cmin cmax]);
cb = colorbar;
cb.Label.String = strrep(fieldName, '_', ' ');
cb.Label.FontSize = 11;

axis equal;
xlim([-0.5 Lx + 0.5]);
ylim([-0.5 Ly + 0.5]);
hold off;
end


%% ========================================================================
% HELPER: Force vector arrows on edges
% =========================================================================
function plotEdgeVectors(celldata, param, edgeForceData, tstep)

Lx = param.Lx;
Ly = param.Ly;
nEdges = length(edgeForceData);

hold on;

% Draw tissue background
plot3Dtissue(celldata.nCells, celldata.r, celldata.connec, param);
rectangle('Position', [0 0 Lx Ly], 'EdgeColor', [0.5 0.5 0.5], 'LineStyle', '--');

% Collect midpoints and force vectors
midX = zeros(nEdges, 1);
midY = zeros(nEdges, 1);
ftX  = zeros(nEdges, 1);
ftY  = zeros(nEdges, 1);
fnX  = zeros(nEdges, 1);
fnY  = zeros(nEdges, 1);

for i = 1:nEdges
    if edgeForceData(i).length < 1e-10
        continue;
    end
    midX(i) = edgeForceData(i).midpoint(1);
    midY(i) = edgeForceData(i).midpoint(2);
    ftX(i)  = edgeForceData(i).force_tangent(1);
    ftY(i)  = edgeForceData(i).force_tangent(2);
    fnX(i)  = edgeForceData(i).force_normal(1);
    fnY(i)  = edgeForceData(i).force_normal(2);
end

% Plot tangential forces (red)
quiver(midX, midY, ftX, ftY, 0.5, 'r', 'LineWidth', 1.2, 'MaxHeadSize', 0.3);

% Plot normal forces (blue)
quiver(midX, midY, fnX, fnY, 0.5, 'b', 'LineWidth', 1.0, 'MaxHeadSize', 0.3);

legend('', 'Tangential (tension)', 'Normal (pressure)', 'Location', 'best');

axis equal;
xlim([-0.5 Lx + 0.5]);
ylim([-0.5 Ly + 0.5]);
hold off;
end


%% ========================================================================
% HELPER: Blue-White-Red colormap
% =========================================================================
function cmap = blueWhiteRed(n)
% Creates a diverging blue-white-red colormap centered at white
if nargin < 1
    n = 256;
end

half = floor(n / 2);
% Blue to white
r1 = linspace(0, 1, half)';
g1 = linspace(0, 1, half)';
b1 = ones(half, 1);

% White to red
r2 = ones(n - half, 1);
g2 = linspace(1, 0, n - half)';
b2 = linspace(1, 0, n - half)';

cmap = [r1 g1 b1; r2 g2 b2];
end


%% ========================================================================
% HELPER: Force vector arrows per unit edge length
% =========================================================================
function plotEdgeVectorsPerLength(celldata, param, edgeForceData, tstep)

Lx = param.Lx;
Ly = param.Ly;
nEdges = length(edgeForceData);

hold on;

% Draw tissue background
plot3Dtissue(celldata.nCells, celldata.r, celldata.connec, param);
rectangle('Position', [0 0 Lx Ly], 'EdgeColor', [0.5 0.5 0.5], 'LineStyle', '--');

% Collect midpoints and per-length force vectors
midX = zeros(nEdges, 1);
midY = zeros(nEdges, 1);
ftX  = zeros(nEdges, 1);
ftY  = zeros(nEdges, 1);
fnX  = zeros(nEdges, 1);
fnY  = zeros(nEdges, 1);

for i = 1:nEdges
    if edgeForceData(i).length < 1e-10
        continue;
    end
    midX(i) = edgeForceData(i).midpoint(1);
    midY(i) = edgeForceData(i).midpoint(2);
    ftX(i)  = edgeForceData(i).force_tangent_per_length(1);
    ftY(i)  = edgeForceData(i).force_tangent_per_length(2);
    fnX(i)  = edgeForceData(i).force_normal_per_length(1);
    fnY(i)  = edgeForceData(i).force_normal_per_length(2);
end

% Compute magnitudes
magT = sqrt(ftX.^2 + ftY.^2);
magN = sqrt(fnX.^2 + fnY.^2);
maxT = max(magT);
maxN = max(magN);

% Manual scaling to avoid quiver auto-scale distortions
% We want the maximum arrow length to be about 0.5 units in the plot
target_length = 0.5;

if maxN > 1e-12
    scaleN = target_length / maxN;
else
    scaleN = 1;
end
fnX_scaled = fnX * scaleN;
fnY_scaled = fnY * scaleN;

if maxT > 1e-12
    scaleT = target_length / maxT;
else
    scaleT = 1;
end
ftX_scaled = ftX * scaleT;
ftY_scaled = ftY * scaleT;

if scaleT > scaleN * 2
    tangLabel = sprintf('Tangential/length (x%.1f)', scaleT / scaleN);
else
    tangLabel = 'Tangential/length';
end

% Plot tangential forces per length (green) - scale = 0 turns OFF auto-scaling
hT = quiver(midX, midY, ftX_scaled, ftY_scaled, 0, 'g', 'LineWidth', 1.5, 'MaxHeadSize', 0.5);

% Plot normal forces per length (blue) - scale = 0 turns OFF auto-scaling
hN = quiver(midX, midY, fnX_scaled, fnY_scaled, 0, 'b', 'LineWidth', 1.2, 'MaxHeadSize', 0.5);

% Use explicit handles so legend colors are correct
legend([hT, hN], {tangLabel, 'Normal/length'}, 'Location', 'best');

axis equal;
xlim([-0.5 Lx + 0.5]);
ylim([-0.5 Ly + 0.5]);
hold off;
end
