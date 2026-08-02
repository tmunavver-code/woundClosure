function [coordinates, connectivity, edgedata,verttocell, T1flagVec,T1relaxstepcountVec, nT1]   = checkT1transitions(nCells,coordinates, connectivity, edgedata,verttocell,param,T1flagVec,T1relaxstepcountVec,imAugmented)
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

if imAugmented == 1
    label = 'augmented';
else
    label = 'reduced';
end

for cellID = 1:nCells
    verticesmat = connectivity{cellID};
    edgelenmat  = edgedata{cellID};
    
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
            if edgelenmat(j) < current_T1_TOL
                
                % Find the surrounding cells to the edge to be flipped
                [cellID1, cellID2, cellID3, cellID4] = getT1cellIDs(vertID1, vertID2,  verttocell,connectivity);
                
                % Perform the T1 transition
                % --- MODIFIED: Pass current_T1_TOL as the last argument ---
                [coordinates, connectivity, edgedata,verttocell] = performT1swap(cellID1, cellID2, cellID3, cellID4, coordinates, connectivity, edgedata,verttocell, vertID1, vertID2, param,imAugmented, current_T1_TOL);
                
                % Display status
                fprintf(1,'-------------\nT1 transition performed in %s \n',label);
                fprintf(1,'V1 = %3d, V2 = %3d\n',vertID1, vertID2);
                fprintf(1,'C1 = %3d, C2 = %3d, C3 = %3d, C4 = %3d]\n',cellID1, cellID2, cellID3, cellID4);
                fprintf(1,'-------------\n');
                
                nT1 = nT1 + 1;
                
                % Set relaxation flag on the VERTICES involved
                T1flagVec(vertID1) = 1;
                T1flagVec(vertID2) = 1;
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