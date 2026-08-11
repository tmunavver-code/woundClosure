function fig = validateTrepatForces(edgeForceData, param, tstep)
% validateTrepatForces: Validates the spatial force field against Trepat 2014.
% Computes radial and azimuthal components of edge forces with respect to
% the wound center, bins them by radial distance, and plots the profile.

Lx = param.Lx;
Ly = param.Ly;
rc = param.woundCircular.center;

nEdges = length(edgeForceData);
radii = zeros(nEdges, 1);
F_rad = zeros(nEdges, 1);
F_azi = zeros(nEdges, 1);

for i = 1:nEdges
    if edgeForceData(i).length < 1e-10
        continue;
    end
    
    mid = edgeForceData(i).midpoint;
    
    % Vector from center to edge midpoint (with periodic wrapping)
    d = mid - rc;
    d(1) = d(1) - Lx * round(d(1) / Lx);
    d(2) = d(2) - Ly * round(d(2) / Ly);
    
    r = norm(d);
    radii(i) = r;
    
    if r > 1e-10
        er = d / r;
        et = [-er(2), er(1)]; % Azimuthal unit vector (CCW)
        
        % Total force vector per unit length on the edge
        F_total = edgeForceData(i).force_tangent_per_length + edgeForceData(i).force_normal_per_length;
        
        % Project onto radial and azimuthal directions
        F_rad(i) = dot(F_total, er);
        F_azi(i) = dot(F_total, et);
    end
end

% Filter out invalid/zero length edges
valid = radii > 1e-10;
r_valid = radii(valid);
Fr_valid = F_rad(valid);
Ft_valid = F_azi(valid);

% Bin the data by radial distance
nBins = 30;
maxR = max(r_valid);
edges = linspace(0, maxR, nBins+1);
binCenters = (edges(1:end-1) + edges(2:end)) / 2;

Fr_binned = zeros(nBins, 1);
Ft_binned = zeros(nBins, 1);

for b = 1:nBins
    inBin = r_valid >= edges(b) & r_valid < edges(b+1);
    if any(inBin)
        Fr_binned(b) = mean(Fr_valid(inBin));
        Ft_binned(b) = mean(Ft_valid(inBin));
    end
end

% Plotting
fig = figure('visible', 'off', 'Position', [100 100 800 500]);
plot(binCenters, Fr_binned, 'r-o', 'LineWidth', 2, 'MarkerFaceColor', 'r');
hold on;
plot(binCenters, Ft_binned, 'b-s', 'LineWidth', 2, 'MarkerFaceColor', 'b');
plot([0 maxR], [0 0], 'k--', 'LineWidth', 1); % Zero line

% Mark wound radius approximately (based on area or known parameter)
% If param.woundCircular.radius exists, plot it
if isfield(param.woundCircular, 'radius')
    xline(param.woundCircular.radius, 'g--', 'LineWidth', 2, 'Label', 'Initial Wound Edge');
end

xlabel('Distance from Wound Center (\mum)');
ylabel('Average Force per Unit Length');
title(sprintf('Spatial Force Field Validation (Trepat 2014) at Timestep %d', tstep));
legend('Radial Force (F_{rad})', 'Azimuthal Force (F_{azi})', 'Location', 'best');
grid on;
hold off;

end
