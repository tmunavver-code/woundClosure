function fig = validateTrepatForces(edgeForceData, param, tstep, varargin)
% VALIDATETREPATFORCES Validates spatial force profiles against Trepat et al.
% Computes radial and azimuthal components of edge forces relative to wound center.
% USAGE: fig = validateTrepatForces(edgeForceData, param, tstep, varargin)

%% -------- Parse options --------
p = inputParser;
p.addParameter('timeAverage', false, @islogical);
p.addParameter('nBins', [], @(x) isempty(x) || (isscalar(x) && x > 0));
p.addParameter('minCount', 3, @(x) isscalar(x) && x >= 1);
p.addParameter('minEdgeFrac', 0.5, @(x) isscalar(x) && x > 0);
p.parse(varargin{:});
opt = p.Results;

Lx = param.Lx;
Ly = param.Ly;
rc = param.woundCircular.center;

% Minimum edge length filter to avoid division by zero near degenerate edges
if isfield(param, 'T1_TOL_bulk')
    minEdgeLen = opt.minEdgeFrac * param.T1_TOL_bulk;
else
    minEdgeLen = 1e-3;  % fallback
end

%% -------- Gather edges (single snapshot or pooled) --------
if opt.timeAverage
    % First arg is a struct array of snapshots, each with a .data field.
    snaps = edgeForceData;
    allEF = [];
    for s = 1:numel(snaps)
        if isfield(snaps(s), 'data')
            allEF = [allEF, snaps(s).data]; %#ok<AGROW>
        else
            allEF = [allEF, snaps(s)];      %#ok<AGROW>  % already an EF array
        end
    end
    edgeForceData = allEF;
end

nEdges = numel(edgeForceData);
radii = nan(nEdges, 1);
F_rad = nan(nEdges, 1);
F_azi = nan(nEdges, 1);

nSkippedShort = 0;

for i = 1:nEdges
    L = edgeForceData(i).length;

    % --- Physical minimum-length guard (Bug 1) ---
    if L < minEdgeLen
        nSkippedShort = nSkippedShort + 1;
        continue;
    end

    mid = edgeForceData(i).midpoint;

    % Vector from wound center to edge midpoint (periodic wrapping)
    d = mid - rc;
    d(1) = d(1) - Lx * round(d(1) / Lx);
    d(2) = d(2) - Ly * round(d(2) / Ly);

    r = norm(d);
    if r <= 1e-10
        continue;  % edge sitting on the center; radial direction undefined
    end
    radii(i) = r;

    er = d / r;
    et = [-er(2), er(1)];   % azimuthal (CCW) unit vector

    % Total force per unit length on this edge
    F_total = edgeForceData(i).force_tangent_per_length + ...
              edgeForceData(i).force_normal_per_length;

    F_rad(i) = dot(F_total, er);
    F_azi(i) = dot(F_total, et);
end

valid   = ~isnan(radii);
r_valid = radii(valid);
Fr_valid = F_rad(valid);
Ft_valid = F_azi(valid);

if isempty(r_valid)
    warning('validateTrepatForces: no valid edges after filtering.');
    fig = figure('visible', 'off');
    return;
end

%% -------- Adaptive bin count (Bug 3) --------
% Aim for ~minCount*2 edges per bin on average; cap between 6 and 15 bins.
if isempty(opt.nBins)
    targetPerBin = max(2 * opt.minCount, 6);
    nBins = floor(numel(r_valid) / targetPerBin);
    nBins = max(6, min(15, nBins));
else
    nBins = opt.nBins;
end

maxR = max(r_valid);
binEdges = linspace(0, maxR, nBins + 1);
binCenters = (binEdges(1:end-1) + binEdges(2:end)) / 2;

Fr_binned  = nan(nBins, 1);
Ft_binned  = nan(nBins, 1);
Fr_err     = nan(nBins, 1);
Ft_err     = nan(nBins, 1);
N_perbin   = zeros(nBins, 1);

for b = 1:nBins
    if b < nBins
        inBin = r_valid >= binEdges(b) & r_valid < binEdges(b+1);
    else
        inBin = r_valid >= binEdges(b) & r_valid <= binEdges(b+1); % include max
    end
    n = sum(inBin);
    N_perbin(b) = n;

    % --- Blank out under-sampled bins as NaN (Bugs 4 & 5) ---
    if n >= opt.minCount
        Fr_binned(b) = mean(Fr_valid(inBin));
        Ft_binned(b) = mean(Ft_valid(inBin));
        % standard error of the mean
        Fr_err(b) = std(Fr_valid(inBin)) / sqrt(n);
        Ft_err(b) = std(Ft_valid(inBin)) / sqrt(n);
    end
    % else leave as NaN so the line breaks instead of drawing phantom zeros
end

%% -------- Plotting --------
fig = figure('visible', 'off', 'Position', [100 100 900 550]);

% Top panel: force profiles with SEM error bars
ax1 = subplot(4, 1, 1:3);
hold(ax1, 'on');
errorbar(binCenters, Fr_binned, Fr_err, 'r-o', 'LineWidth', 2, ...
    'MarkerFaceColor', 'r', 'CapSize', 4, 'DisplayName', 'Radial Force (F_{rad})');
errorbar(binCenters, Ft_binned, Ft_err, 'b-s', 'LineWidth', 2, ...
    'MarkerFaceColor', 'b', 'CapSize', 4, 'DisplayName', 'Azimuthal Force (F_{azi})');
plot([0 maxR], [0 0], 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');

if isfield(param.woundCircular, 'radius')
    xline(param.woundCircular.radius, 'g--', 'LineWidth', 2, ...
        'Label', 'Initial Wound Edge', 'HandleVisibility', 'off');
end

ylabel('Avg Force per Unit Length');
if opt.timeAverage
    ttl = sprintf('Spatial Force Field vs Trepat 2014 (time-averaged, %d bins)', nBins);
else
    ttl = sprintf('Spatial Force Field vs Trepat 2014 (t=%d, %d bins)', tstep, nBins);
end
title(ttl);
legend('Location', 'best');
grid(ax1, 'on');
xlim([0 maxR]);
hold(ax1, 'off');

% Bottom panel: edges per bin, so you can see where the signal is trustworthy
ax2 = subplot(4, 1, 4);
bar(binCenters, N_perbin, 'FaceColor', [0.6 0.6 0.6], 'EdgeColor', 'none');
hold(ax2, 'on');
yline(opt.minCount, 'r--', 'LineWidth', 1, ...
    'Label', sprintf('min=%d', opt.minCount));
if isfield(param.woundCircular, 'radius')
    xline(param.woundCircular.radius, 'g--', 'LineWidth', 2, 'HandleVisibility', 'off');
end
ylabel('# edges');
xlabel('Distance from Wound Center (\mum)');
grid(ax2, 'on');
xlim([0 maxR]);
hold(ax2, 'off');

% Console diagnostics
fprintf('validateTrepatForces: %d valid edges, %d bins, %d short edges excluded (< %.4f).\n', ...
    numel(r_valid), nBins, nSkippedShort, minEdgeLen);
nBlank = sum(N_perbin < opt.minCount);
if nBlank > 0
    fprintf('  %d/%d bins blanked (fewer than %d edges) -- tissue is under-sampled.\n', ...
        nBlank, nBins, opt.minCount);
    if ~opt.timeAverage
        fprintf('  Consider re-running with ''timeAverage'', true over your snapshots.\n');
    end
end
end
