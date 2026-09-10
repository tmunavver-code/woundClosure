function [celldata, T1flagVec, intercalation_happened, param] = checkWoundIntercalations(celldata, param, T1flagVec, edgeForceData, tstep)
%
% VERSION 2.1 (Corrected):
% - This is the "complete detachment" logic.
% - It now returns the modified 'param' struct so the main
%   script can update the 'cellIDtoContract' list for plotting.
%

% Flag to tell the main loop if we need to recalculate geometry
intercalation_happened = false;

% Get data from celldata struct
coordinates = celldata.r;
connectivity = celldata.connec;
verttocell = celldata.verttocell;
woundEdges = celldata.woundEdges;
woundCells = param.cellIDtoContract; % List of all wound cells

% Energy BEFORE any change this call, used as the gate below. Computed
% once since only one intercalation is now committed per call (see the
% one-per-call fix below), so there is no risk of it going stale mid-loop.
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

        % --- 2. IDENTIFY ALL PLAYERS ---
        % REDESIGNED (matching Tetley et al. 2019's own stated method,
        % confirmed against their Supplementary Video 4: "the wound itself
        % was treated as a cell, meaning an intercalating tetrad was
        % formed of three wound-edge cells and the wound"). This is now a
        % genuine T1-style neighbor swap with the wound playing the role
        % of the 4th "cell" -- NOT the old merge-and-strip approach, which
        % progressively shrank the margin cell toward degeneracy. In the
        % video, a cell pushed off the margin remains a fully intact,
        % normal cell elsewhere in the tissue; it never shrinks or
        % vanishes from a single event.
        %
        %   cellA: the real cell bordering BOTH v1_master and v2_orphan
        %          (touches the wound along this edge). Loses this edge --
        %          pushed off the margin -- but stays a complete cell.
        %   cellB: the other real cell at v1_master (not v2_orphan).
        %          Gains v2_orphan.
        %   cellC: the other real cell at v2_orphan (not v1_master).
        %          Gains v1_master.
        cells_at_v1 = verttocell{v1_master};
        cells_at_v2 = verttocell{v2_orphan};
        cellA_candidates = intersect(cells_at_v1, cells_at_v2);
        cellA = intersect(cellA_candidates, woundCells);

        if length(cellA) ~= 1
            continue;
        end

        % cellA must have at least 4 vertices -- it is about to lose one
        % (v2_orphan) and a real cell must never drop below 3.
        if length(connectivity{cellA}) < 4
            continue;
        end

        % cellB/cellC can legitimately not exist: if a vertex's OTHER edge
        % (besides this candidate wound edge) is ALSO a wound edge, cellA
        % is the only real cell touching it -- a genuine "corner"/spike of
        % the margin poking into the wound. A 4-cell-style handoff has no
        % partner to give a vertex to there. This is exactly what was
        % freezing margin turnover once the shrinking wound made most
        % remaining vertices corners rather than smooth 2-neighbor points
        % (empirically: margin stuck at 13 cells from step ~300 onward in
        % a 1000-step run, even as wound area kept falling toward zero).
        %
        % Fix: handle it as a second, distinct topological case -- trim
        % the corner vertex out of cellA entirely (no handoff needed,
        % since both its neighboring faces are already the wound; this
        % merges the two wound edges meeting there into one). Still
        % gated by the same force-based Bell's-law trigger above and the
        % energy/Metropolis gate below -- only the geometric mechanics
        % differ between the two cases.
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
            % Both ends are corners, or genuinely ambiguous topology --
            % skip rather than guess; a neighboring edge will usually
            % resolve this from the other side instead.
            continue;
        end

        % --- 3. PERFORM THE SWAP (trial) ---
        % Trial on local copies, gated on energy exactly as before -- only
        % the mechanics of what a successful flip does have changed.
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
            % Metropolis-style acceptance (see checkT1transitions.m for
            % the same treatment and rationale) -- occasional uphill
            % detachments allowed instead of a strict always-downhill rule.
            P_accept_uphill = exp(-dE / param.T_eff_T2);
            accept_intercalation = (rand() < P_accept_uphill);
        else
            accept_intercalation = false;
        end

        if ~accept_intercalation
            % Rejected: leaves coordinates/connectivity/verttocell untouched,
            % try the next candidate wound edge instead.
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

        % BUG FIX: process at most one intercalation per call. The
        % original loop could commit many detachments per timestep off of
        % a single stale copy of woundEdges/connectivity captured at the
        % top of the function -- the same staleness issue
        % README_force_based_transitions.md already flags for
        % checkT1transitions.m ("Phase 1b", stale forces within a sweep).
        % Empirically this is what let the margin grow instead of shrink.
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