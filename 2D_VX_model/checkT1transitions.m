function [coordinates, connectivity, edgedata,verttocell, T1flagVec,T1relaxstepcountVec, nT1, T1ForceLog]   = checkT1transitions(nCells,coordinates, connectivity, edgedata,verttocell,param,T1flagVec,T1relaxstepcountVec,imAugmented, edgeForceData, tstep, T1ForceLog, celldata)
%% checkT1transitions: Checks for and performs a T1 transition
%
% MODIFIED:
% 1. Implements enhanced T1 intercalation (fluidity) at the wound margin
%    by using a separate, more permissive tolerance (param.T1_TOL_wound).
% 2. Fixes bug in T1flagVec logic: flags are now correctly applied to
%    vertices (vertID1, vertID2) instead of the local edge index (j).
% 3. Passes the triggered T1 tolerance ('current_T1_TOL') into
%    performT1swap.
%

nT1 = param.nT1;

if nargin < 12
    T1ForceLog = [];
end

if imAugmented == 1
    label = 'augmented';
else
    label = 'reduced';
end

T1_happened_this_step = false;

for cellID = 1:nCells
    if T1_happened_this_step
        break;
    end
    
    verticesmat = connectivity{cellID};
    edgelenmat  = edgedata{cellID};
    
    % --- Guard: skip cells where edgedata is out of sync with vertices ---
    % This can happen after wound intercalations modify connectivity mid-step
    if length(verticesmat) < 3 || length(edgelenmat) < length(verticesmat)
        continue;
    end
    
    for j = 1:length(verticesmat)
        
        % --- Get vertex IDs for the current edge FIRST ---
        vertID1 = verticesmat(j);
        if j < length(verticesmat)
            vertID2 = verticesmat(j+1);
        else % Final edge
            vertID2 = verticesmat(1);
        end
        
        % --- Check T1 relaxation flag on VERTICES (Bug Fix) ---
        % Only proceed if neither vertex is in a relaxation period
        if T1flagVec(vertID1) == 0 && T1flagVec(vertID2) == 0
            
            % --- NEW: Determine correct T1 tolerance ---
            % Start with the default bulk tolerance
            current_T1_TOL = param.T1_TOL_bulk;
            
            % If enhancement is on, check if this cell is a wound cell
            if param.enhanceT1atWound
                % param.cellIDtoContract holds the list of wound edge cells
                if ismember(cellID, param.cellIDtoContract)
                    % Use the more permissive tolerance for wound-edge cells
                    current_T1_TOL = param.T1_TOL_wound;
                end
            end
            % --- END NEW ---
            
            
            % --- Now, check edge length with the correct tolerance ---
            % AND check Force-Based criterion if enabled
            do_flip = false;

            % BUG FIX: previously, when useForceBasedTransitions was true but
            % the caller didn't pass edgeForceData (nargin<10, e.g. the
            % myMainClassical.m call sites, which only pass 9 args), do_flip
            % could never become true -- T1 was silently dead code. Now falls
            % back to the kinematic rule whenever force data isn't available,
            % matching the fallback checkWoundIntercalations.m already has.
            if isfield(param, 'useForceBasedTransitions') && param.useForceBasedTransitions && nargin >= 10 && ~isempty(edgeForceData)
                % Force-Based logic: length is still a NECESSARY pre-condition.
                % Force modulates the probability, but the edge must be short first.
                v_min = min(vertID1, vertID2);
                v_max = max(vertID1, vertID2);
                idx = find([edgeForceData.v1] == v_min & [edgeForceData.v2] == v_max, 1);
                if ~isempty(idx)
                    edge_len = edgeForceData(idx).length;

                    % GATE 1: Edge must be shorter than the tolerance (length pre-condition)
                    if edge_len < current_T1_TOL
                        force_tangent_per_length = edgeForceData(idx).total_tension / edge_len;

                        % Length/crowding modulation: f_beta decreases as edge shortens,
                        % causing the probability to self-accelerate near collapse.
                        % T1 gates on TANGENTIAL force -- uses the T1-specific scale.
                        f_beta_eff = max(param.f_beta_T1 - param.f_beta_length_slope_T1 * (current_T1_TOL - edge_len), 0.01);

                        % GATE 2: Bell's-law stochastic probability
                        F_drive = force_tangent_per_length;
                        k_off = param.k_off0 * exp(F_drive / f_beta_eff);
                        P_flip = 1 - exp(-param.deltat * k_off);

                        if rand() < P_flip
                            do_flip = true;
                        end
                    end
                end
            else
                % Kinematic (Length-Based) logic
                if edgelenmat(j) < current_T1_TOL
                    do_flip = true;
                end
            end
            
            if do_flip
                
                % --- Log forces before flip ---
                v_min = min(vertID1, vertID2);
                v_max = max(vertID1, vertID2);
                if nargin >= 10 && ~isempty(edgeForceData)
                    idx = find([edgeForceData.v1] == v_min & [edgeForceData.v2] == v_max, 1);
                    if ~isempty(idx)
                        tension_per_L = edgeForceData(idx).total_tension / edgeForceData(idx).length;
                        pressure = edgeForceData(idx).pressure_area;
                        % log: [tstep, v1, v2, length, tension_per_length, pressure]
                        if nargin < 12 || isempty(T1ForceLog)
                            T1ForceLog = [tstep, v_min, v_max, edgelenmat(j), tension_per_L, pressure];
                        else
                            T1ForceLog(end+1, :) = [tstep, v_min, v_max, edgelenmat(j), tension_per_L, pressure];
                        end
                    end
                end
                
                % --- Energy Gate pre-calculation ---
                if nargin >= 13 && ~isempty(celldata)
                    E_old = getTissueEnergyClassical(celldata, param);
                else
                    E_old = 0; % Fallback if celldata isn't provided
                end
                
                % Find the surrounding cells to the edge to be flipped
                [cellID1, cellID2, cellID3, cellID4] = getT1cellIDs(vertID1, vertID2,  verttocell,connectivity);
                
                % Guard: skip if we couldn't find 4 valid surrounding cells
                if cellID1 == -1 || cellID2 == -1 || cellID3 == -1 || cellID4 == -1
                    continue;
                end
                
                % Backup state before swap
                coords_backup = coordinates;
                conn_backup = connectivity;
                edge_backup = edgedata;
                vtc_backup = verttocell;
                
                % Perform the T1 transition
                % --- MODIFIED: Pass current_T1_TOL as the last argument ---
                [coordinates, connectivity, edgedata,verttocell] = performT1swap(cellID1, cellID2, cellID3, cellID4, coordinates, connectivity, edgedata,verttocell, vertID1, vertID2, param,imAugmented, current_T1_TOL);
                
                % --- Energy Gate post-calculation ---
                accept_swap = true;
                if nargin >= 13 && ~isempty(celldata)
                    celldata_temp = celldata;
                    celldata_temp.r = coordinates;
                    celldata_temp.connec = connectivity;
                    celldata_temp.EdgeData = edgedata;
                    celldata_temp.verttocell = verttocell;
                    
                    % Need to update Areas and Perimeters for the new geometry
                    celldata_temp.A = getCellAreas(celldata_temp, param);
                    [celldata_temp.P, celldata_temp.EdgeData] = getCellPerimeters(celldata_temp.nCells, celldata_temp.r, celldata_temp.connec, param, 0);
                    
                    E_new = getTissueEnergyClassical(celldata_temp, param);
                    dE = E_new - E_old;

                    if dE < 0
                        accept_swap = true;
                    elseif isfield(param, 'T_eff_T1') && param.T_eff_T1 > 0
                        % Metropolis-style acceptance: allow occasional
                        % uphill moves instead of a strict always-downhill
                        % rule, so the system can escape a locally-stuck
                        % configuration via fluctuation and relax down
                        % afterward (the way real active/thermal systems,
                        % and Kaurin-Arroyo's stochastic bond-breaking,
                        % can transiently do). T_eff_T1 sets the scale of
                        % uphill jump this system tolerates -- calibrated
                        % from this model's own logged dE, same approach
                        % as f_beta_T1/f_beta_T2.
                        P_accept_uphill = exp(-dE / param.T_eff_T1);
                        accept_swap = (rand() < P_accept_uphill);
                    else
                        % No T_eff_T1 configured -- fall back to the old
                        % strict always-downhill rule.
                        accept_swap = false;
                    end
                end
                
                if accept_swap
                    % Display status
                    fprintf(1,'-------------\nT1 transition performed in %s \n',label);
                    fprintf(1,'V1 = %3d, V2 = %3d\n',vertID1, vertID2);
                    fprintf(1,'C1 = %3d, C2 = %3d, C3 = %3d, C4 = %3d]\n',cellID1, cellID2, cellID3, cellID4);
                    if nargin >= 13 && ~isempty(celldata)
                        fprintf(1,'Energy change: %f\n', E_new - E_old);
                    end
                    fprintf(1,'-------------\n');
                    
                    nT1 = nT1 + 1;
                    
                    % Set relaxation flag on the VERTICES involved
                    T1flagVec(vertID1) = 1;
                    T1flagVec(vertID2) = 1;
                    
                    T1_happened_this_step = true;
                    break; % Break out of inner loop over edges
                else
                    % Revert the swap
                    if nargin >= 13 && ~isempty(celldata)
                        fprintf(1, 'T1_REJECTED_DE: %f\n', dE);
                    end
                    coordinates = coords_backup;
                    connectivity = conn_backup;
                    edgedata = edge_backup;
                    verttocell = vtc_backup;
                end
            end
            
        else
            % Do nothing, vertices are in relaxation period
        end
    end
end

% --- This relaxation loop is correct ---
% It iterates over the global vertex list and counts down
% any flagged vertices.
for i = 1:length(T1flagVec)
    if T1flagVec(i) == 1
        T1relaxstepcountVec(i) = T1relaxstepcountVec(i) + 1;
    end
    
    if T1relaxstepcountVec(i) > param.maxT1relaxsteps
        T1flagVec(i) = 0;
        T1relaxstepcountVec(i) = 0;
    end
end

end