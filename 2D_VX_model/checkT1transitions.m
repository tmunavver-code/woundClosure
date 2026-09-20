function [coordinates, connectivity, edgedata,verttocell, T1flagVec,T1relaxstepcountVec, nT1, T1ForceLog]   = checkT1transitions(nCells,coordinates, connectivity, edgedata,verttocell,param,T1flagVec,T1relaxstepcountVec,imAugmented, edgeForceData, tstep, T1ForceLog, celldata)
%% checkT1transitions: Evaluates and executes T1 neighbor exchanges
% Uses Bell's law on edge tension and Metropolis uphill energy acceptance.

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
    
    % Guard: skip cells where edgedata is out of sync with vertices
    if length(verticesmat) < 3 || length(edgelenmat) < length(verticesmat)
        continue;
    end
    
    for j = 1:length(verticesmat)
        
        vertID1 = verticesmat(j);
        if j < length(verticesmat)
            vertID2 = verticesmat(j+1);
        else
            vertID2 = verticesmat(1);
        end
        
        % Check relaxation flag on vertices
        if T1flagVec(vertID1) == 0 && T1flagVec(vertID2) == 0
            
            % Determine T1 tolerance (bulk vs wound margin)
            current_T1_TOL = param.T1_TOL_bulk;
            if param.enhanceT1atWound
                if ismember(cellID, param.cellIDtoContract)
                    current_T1_TOL = param.T1_TOL_wound;
                end
            end
            
            do_flip = false;

            % Force-based or kinematic logic
            if isfield(param, 'useForceBasedTransitions') && param.useForceBasedTransitions && nargin >= 10 && ~isempty(edgeForceData)
                v_min = min(vertID1, vertID2);
                v_max = max(vertID1, vertID2);
                idx = find([edgeForceData.v1] == v_min & [edgeForceData.v2] == v_max, 1);
                if ~isempty(idx)
                    edge_len = edgeForceData(idx).length;

                    % Edge must be shorter than tolerance
                    if edge_len < current_T1_TOL
                        force_tangent_per_length = edgeForceData(idx).total_tension / edge_len;

                        % Modulate f_beta with edge shortening
                        f_beta_eff = max(param.f_beta_T1 - param.f_beta_length_slope_T1 * (current_T1_TOL - edge_len), 0.01);

                        % Bell's law stochastic probability
                        F_drive = force_tangent_per_length;
                        k_off = param.k_off0 * exp(F_drive / f_beta_eff);
                        P_flip = 1 - exp(-param.deltat * k_off);

                        if rand() < P_flip
                            do_flip = true;
                        end
                    end
                end
            else
                % Kinematic (length-based) logic
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
                        % Metropolis uphill acceptance
                        P_accept_uphill = exp(-dE / param.T_eff_T1);
                        accept_swap = (rand() < P_accept_uphill);
                    else
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
                    
                    % Set relaxation flag on vertices
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
            % Vertices in relaxation period
        end
    end
end

% Decrement relaxation counter on flagged vertices
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