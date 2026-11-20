function [celldata, param, seal_happened, eventInfo] = checkWoundSealing(celldata, param, tstep)
% CHECKWOUNDSEALING Merges opposing wound margin vertices upon contact.
% Gated by distance threshold, Bell's law closing force, and Metropolis energy.

seal_happened = false;
eventInfo = struct('timestep', tstep, 'happened', false, 'reason', '', ...
    'v1', NaN, 'v2', NaN, 'separation', NaN, 'forceDrive', NaN, ...
    'bellProbability', NaN, 'energyChange', NaN, ...
    'boundaryCellsBefore', NaN, 'boundaryCellsAfter', NaN, ...
    'mergedVertexDegree', NaN);

if ~isfield(param, 'enableWoundSealing') || ~param.enableWoundSealing
    return;
end
if ~isfield(celldata, 'woundEdges') || isempty(celldata.woundEdges)
    eventInfo.reason = 'no_wound_edges';
    return;
end

boundaryCellsBefore = numel(unique(param.cellIDtoContract));
eventInfo.boundaryCellsBefore = boundaryCellsBefore;
if isfield(param, 'maxBoundaryCellsForSealing') && ...
        boundaryCellsBefore > param.maxBoundaryCellsForSealing
    eventInfo.reason = 'waiting_for_boundary_handoffs';
    return;
end

Lx = param.Lx;
Ly = param.Ly;

marginVerts = unique(celldata.woundEdges(:));
if numel(marginVerts) < 4
    return;
end

E_old = getTissueEnergyClassical(celldata, param);
if isfield(param, 'bellProbabilityThresholdSeal')
    probabilityThreshold = param.bellProbabilityThresholdSeal;
else
    probabilityThreshold = param.bellProbabilityThreshold;
end

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

        if P_seal < probabilityThreshold
            continue;
        end

        % --- Trial merge ---
        [coords_trial, connec_trial, ~, ok] = performWoundSeal( ...
            v1, v2, celldata.r, celldata.connec, celldata.verttocell);
        if ~ok
            continue;
        end

        % --- Gate 3: energy ---
        celldata_temp = celldata;
        celldata_temp.r = coords_trial;
        celldata_temp.connec = connec_trial;
        celldata_temp.verttocell = getVertexToCellIDMap(size(coords_trial,1), ...
            celldata_temp.nCells, connec_trial);
        celldata_temp.A = getCellAreas(celldata_temp, param);
        [celldata_temp.P, celldata_temp.EdgeData] = getCellPerimeters(celldata_temp.nCells, celldata_temp.r, celldata_temp.connec, param, 0);
        [newWoundCells, newWoundEdges] = findWoundEdgeCells(celldata_temp);
        celldata_temp.woundEdges = newWoundEdges;
        celldata_temp.cellIDtoContract = newWoundCells;
        param_temp = param;
        param_temp.cellIDtoContract = newWoundCells;

        [topologyValid, ~, maxDegree] = validateTissueTopology(celldata_temp);
        [~, loopValid] = getOrderedWoundLoop(newWoundEdges);
        if isempty(newWoundEdges)
            loopValid = true;
        end
        mergedDegree = numel(celldata_temp.verttocell{v1});
        if ~topologyValid || ~loopValid || ...
                mergedDegree > param.maxSealVertexDegree || ...
                maxDegree > max(param.maxSealVertexDegree, ...
                max(cellfun(@(x) numel(unique(x)), celldata.verttocell))) || ...
                numel(newWoundCells) > boundaryCellsBefore
            eventInfo.reason = 'topology_guard';
            continue;
        end

        E_new = getTissueEnergyClassical(celldata_temp, param_temp);
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
        param = param_temp;
        seal_happened = true;
        param.nSeals = param.nSeals + 1;
        eventInfo.happened = true;
        eventInfo.reason = 'accepted';
        eventInfo.v1 = v1;
        eventInfo.v2 = v2;
        eventInfo.separation = sep;
        eventInfo.forceDrive = F_drive;
        eventInfo.bellProbability = P_seal;
        eventInfo.energyChange = dE;
        eventInfo.boundaryCellsAfter = numel(newWoundCells);
        eventInfo.mergedVertexDegree = mergedDegree;
        fprintf(1, '--- Wound Sealed --- vertices %d and %d merged (gap %.4f, F_drive %.4f)\n', v1, v2, sep, F_drive);
        fprintf(1, 'Energy change: %f\n', dE);

        % One seal per call to avoid stale geometry snapshot
        return;
    end
end

end
