function [celldata, param, seal_happened] = checkWoundSealing(celldata, param, tstep)
%CHECKWOUNDSEALING Seal the wound where opposing margins have come into
%contact -- the final step of closure, which retires cells from the wound
%boundary.
%
% Gated exactly like the other three topological event types in this
% model: a distance pre-condition, then a force-based Bell's-law
% stochastic draw, then a Metropolis energy gate.
%
%   Gate 1 (geometry): two wound-margin vertices closer than
%       param.L_seal_thresh that do NOT already share a cell. Requiring
%       no shared cell is what makes this a merge ACROSS the wound gap
%       (two different sides coming together) rather than a cell pinching
%       itself in half, and it also guarantees no cell ends up listing
%       the surviving vertex twice.
%
%   Gate 2 (force, Bell's law): F_drive is the closing force -- the
%       component of each vertex's own net force directed at the other.
%       Positive means the mechanics are actively pressing the two
%       margins together, which is exactly when adhesion should form.
%       Uses celldata.f, the same force field that moves the vertices, so
%       the gate cannot disagree with the dynamics.
%
%   Gate 3 (energy): accept if the merge lowers tissue energy, otherwise
%       accept with probability exp(-dE/T_eff_seal), as for T1/T2/division.

seal_happened = false;

if ~isfield(param, 'enableWoundSealing') || ~param.enableWoundSealing
    return;
end
if ~isfield(celldata, 'woundEdges') || isempty(celldata.woundEdges)
    return;
end

Lx = param.Lx;
Ly = param.Ly;

marginVerts = unique(celldata.woundEdges(:));
if numel(marginVerts) < 4
    return;
end

E_old = getTissueEnergyClassical(celldata, param);

for a = 1:numel(marginVerts)-1
    v1 = marginVerts(a);
    if isempty(celldata.verttocell{v1}), continue; end

    for b = a+1:numel(marginVerts)
        v2 = marginVerts(b);
        if isempty(celldata.verttocell{v2}), continue; end

        % --- Gate 1a: must not already share a cell (see header) ---
        if ~isempty(intersect(celldata.verttocell{v1}, celldata.verttocell{v2}))
            continue;
        end

        % --- Gate 1b: separation below threshold (periodic) ---
        d = celldata.r(v2,:) - celldata.r(v1,:);
        d(1) = d(1) - Lx * round(d(1) / Lx);
        d(2) = d(2) - Ly * round(d(2) / Ly);
        sep = norm(d);
        if sep >= param.L_seal_thresh || sep < 1e-12
            continue;
        end

        % --- Gate 2: Bell's law on the closing force ---
        u = d / sep;  % unit vector from v1 toward v2
        F_drive = dot(celldata.f(v1,:), u) + dot(celldata.f(v2,:), -u);

        % Self-accelerating as the gap shrinks, same form as T1/T2.
        f_beta_eff = max(param.f_beta_seal - param.f_beta_seal_slope * (param.L_seal_thresh - sep), 0.01);

        k_seal = param.k_seal0 * exp(F_drive / f_beta_eff);
        P_seal = 1 - exp(-param.deltat * k_seal);

        if rand() >= P_seal
            continue;
        end

        % --- Trial merge ---
        [coords_trial, connec_trial, vtc_trial, ok] = performWoundSeal( ...
            v1, v2, celldata.r, celldata.connec, celldata.verttocell);
        if ~ok
            continue;
        end

        % --- Gate 3: energy ---
        celldata_temp = celldata;
        celldata_temp.r = coords_trial;
        celldata_temp.connec = connec_trial;
        celldata_temp.verttocell = vtc_trial;
        celldata_temp.A = getCellAreas(celldata_temp, param);
        [celldata_temp.P, celldata_temp.EdgeData] = getCellPerimeters(celldata_temp.nCells, celldata_temp.r, celldata_temp.connec, param, 0);
        E_new = getTissueEnergyClassical(celldata_temp, param);
        dE = E_new - E_old;

        if dE < 0
            accept = true;
        elseif isfield(param, 'T_eff_seal') && param.T_eff_seal > 0
            accept = (rand() < exp(-dE / param.T_eff_seal));
        else
            accept = false;
        end

        if ~accept
            fprintf(1, 'SEAL_REJECTED_DE: %f\n', dE);
            continue;
        end

        celldata = celldata_temp;
        seal_happened = true;
        param.nSeals = param.nSeals + 1;
        fprintf(1, '--- Wound Sealed --- vertices %d and %d merged (gap %.4f, F_drive %.4f)\n', v1, v2, sep, F_drive);
        fprintf(1, 'Energy change: %f\n', dE);

        % One per call, same rule as the other event types -- avoids
        % compounding several merges off one stale geometry snapshot.
        [celldata, param] = refreshWoundMargin(celldata, param);
        return;
    end
end

end
