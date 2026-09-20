function [celldata, T1flagVec, intercalation_happened, param] = checkWoundIntercalations(celldata, param, T1flagVec, edgeForceData, tstep)
% CHECKWOUNDINTERCALATIONS Evaluates cell detachment and intercalations at the wound
% margin using Bell's law on normal pressure and Metropolis energy gating.

% Flag to tell the main loop if we need to recalculate geometry
intercalation_happened = false;

% Get data from celldata struct
coordinates = celldata.r;
connectivity = celldata.connec;
verttocell = celldata.verttocell;
woundEdges = celldata.woundEdges;
woundCells = param.cellIDtoContract; % List of all wound cells

% Reference tissue energy before candidate event
E_old = getTissueEnergyClassical(celldata, param);

edges_to_remove_idx = []; % To track which wound edges get merged

for i = 1:size(woundEdges, 1)
    
    v1_master = woundEdges(i, 1);
    v2_orphan = woundEdges(i, 2);
    
    % Check if this edge is already being relaxed from a T1 swap
    if T1flagVec(v1_master) == 1 || T1flagVec(v2_orphan) == 1
        continue;
    end
    
    % --- 1. CHECK TRIGGER ---
    vec = coordinates(v2_orphan,:) - coordinates(v1_master,:);
    vec(1) = vec(1) - param.Lx * round(vec(1) / param.Lx);
    vec(2) = vec(2) - param.Ly * round(vec(2) / param.Ly);
    L = norm(vec);
    
    do_flip = false;
    
    if isfield(param, 'useForceBasedTransitions') && param.useForceBasedTransitions
        % Force-Based logic: length is still a NECESSARY pre-condition.
        % The wound edge must already be short — force raises/lowers the probability.
        if L < param.L_intercalation_thresh
            v_min = min(v1_master, v2_orphan);
            v_max = max(v1_master, v2_orphan);
            
            if nargin >= 4 && ~isempty(edgeForceData)
                idx = find([edgeForceData.v1] == v_min & [edgeForceData.v2] == v_max, 1);
                if ~isempty(idx)
                    force_normal_per_length = edgeForceData(idx).pressure_area;
                    
                    % F_drive = normal pressure per length (peeling force)
                    F_drive = force_normal_per_length;
                    
                    % f_beta decreases as L shrinks → self-accelerating near collapse.
                    % T2 gates on NORMAL force (peeling) -- uses the T2-specific scale.
                    f_beta_eff = max(param.f_beta_T2 - param.f_beta_length_slope_T2 * (param.L_intercalation_thresh - L), 0.01);
                    
                    k_off = param.k_off0 * exp(F_drive / f_beta_eff);
                    P_flip = 1 - exp(-param.deltat * k_off);
                    
                    if rand() < P_flip
                        do_flip = true;
                    end
                else
                    % Edge not found in force data — fall back to kinematic rule
                    do_flip = true;
                end
            else
                do_flip = true;  % No force data available — use kinematic rule
            end
        end
    else
        % Kinematic logic
        if L < param.L_intercalation_thresh
            do_flip = true;
        end
    end
    
    if do_flip

        % --- 2. Identify interacting cells around margin edge ---
        % cellA borders both vertices; cellB and cellC border v1 and v2 respectively.
        cells_at_v1 = verttocell{v1_master};
        cells_at_v2 = verttocell{v2_orphan};
        cellA_candidates = intersect(cells_at_v1, cells_at_v2);
        cellA = intersect(cellA_candidates, woundCells);

        if length(cellA) ~= 1
            continue;
        end

        % cellA must have at least 4 vertices to lose one and stay valid
        if length(connectivity{cellA}) < 4
            continue;
        end

        % Check whether handoff or corner trimming applies
        cellB_candidates = setdiff(cells_at_v1, cellA);
        cellC_candidates = setdiff(cells_at_v2, cellA);

        if length(cellB_candidates) == 1 && length(cellC_candidates) == 1 && cellB_candidates(1) ~= cellC_candidates(1)
            mode = 'handoff';
            cellB = cellB_candidates(1);
            cellC = cellC_candidates(1);
        elseif isempty(cellB_candidates) && length(cellC_candidates) == 1
            mode = 'trim_v1'; % v1_master is the corner
        elseif isempty(cellC_candidates) && length(cellB_candidates) == 1
            mode = 'trim_v2'; % v2_orphan is the corner
        else
            % Both ends are corners or ambiguous; skip
            continue;
        end

        % --- 3. Perform trial swap or corner trim ---
        coords_trial = coordinates;
        connec_trial = connectivity;
        vtc_trial = verttocell;

        switch mode
            case 'handoff'
                [coords_trial, connec_trial, vtc_trial] = performWoundMarginT1( ...
                    cellA, cellB, cellC, coords_trial, connec_trial, vtc_trial, ...
                    v1_master, v2_orphan, param);
            case 'trim_v1'
                connec_trial = trimWoundCorner(cellA, connec_trial, v1_master);
                vtc_trial{v1_master} = [];
            case 'trim_v2'
                connec_trial = trimWoundCorner(cellA, connec_trial, v2_orphan);
                vtc_trial{v2_orphan} = [];
        end

        % --- Energy gate ---
        celldata_temp = celldata;
        celldata_temp.r = coords_trial;
        celldata_temp.connec = connec_trial;
        celldata_temp.verttocell = vtc_trial;
        celldata_temp.A = getCellAreas(celldata_temp, param);
        [celldata_temp.P, celldata_temp.EdgeData] = getCellPerimeters(celldata_temp.nCells, celldata_temp.r, celldata_temp.connec, param, 0);
        E_new = getTissueEnergyClassical(celldata_temp, param);
        dE = E_new - E_old;

        if dE < 0
            accept_intercalation = true;
        elseif isfield(param, 'T_eff_T2') && param.T_eff_T2 > 0
            % Metropolis uphill acceptance
            P_accept_uphill = exp(-dE / param.T_eff_T2);
            accept_intercalation = (rand() < P_accept_uphill);
        else
            accept_intercalation = false;
        end

        if ~accept_intercalation
            % Rejected trial
            fprintf(1, 'T2_REJECTED_DE: %f\n', dE);
            continue;
        end

        % --- Commit the trial ---
        coordinates = coords_trial;
        connectivity = connec_trial;
        verttocell = vtc_trial;
        intercalation_happened = true;

        % E. Mark Edge for Removal
        edges_to_remove_idx = [edges_to_remove_idx, i]; %#ok<AGROW>

        % F. Set Relaxation Flag
        T1flagVec(v1_master) = 1;

        switch mode
            case 'handoff'
                fprintf(1,'--- Wound-Margin T1 Intercalation Performed ---\n');
                fprintf(1,'Cell %d pushed off the margin (wound-edge %d-%d), neighbors %d and %d take over.\n', cellA, v1_master, v2_orphan, cellB, cellC);
            case 'trim_v1'
                fprintf(1,'--- Wound-Margin Corner Trimmed ---\n');
                fprintf(1,'Vertex %d (corner of cell %d) absorbed into the wound.\n', v1_master, cellA);
            case 'trim_v2'
                fprintf(1,'--- Wound-Margin Corner Trimmed ---\n');
                fprintf(1,'Vertex %d (corner of cell %d) absorbed into the wound.\n', v2_orphan, cellA);
        end
        fprintf(1,'Energy change: %f\n', E_new - E_old);

        % Process at most one intercalation per call to avoid stale geometry compounding
        break;
    end
end

% --- 4. CLEAN UP & UPDATE STATE ---
if intercalation_happened
    % Commit connectivity & coordinate changes
    celldata.r = coordinates;
    celldata.connec = connectivity;
    celldata.verttocell = verttocell;
    
    % The wound margin is now different. (THIS IS THE SLOW STEP)
    [newWoundCells, newWoundEdges] = findWoundEdgeCells(celldata);
    
    celldata.woundEdges = newWoundEdges;
    param.cellIDtoContract = newWoundCells; % Update list for plotting
    celldata.cellIDtoContract = newWoundCells; % Also store in celldata
end

end