function [celldata, T1flagVec, intercalation_happened, param] = checkWoundIntercalations(celldata, param, T1flagVec)
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
    
    if L < param.L_intercalation_thresh
        
        % --- 2. IDENTIFY ALL PLAYERS ---
        cells_at_v1 = verttocell{v1_master};
        cells_at_v2 = verttocell{v2_orphan};
        Cell_W_candidates = intersect(cells_at_v1, cells_at_v2);
        Cell_W = intersect(Cell_W_candidates, woundCells);
        
        if length(Cell_W) ~= 1
            continue;
        end
        
        cells_to_rewire = setdiff(cells_at_v2, cells_at_v1);
        
        % --- 3. PERFORM VERTEX MERGE & ORPHAN ---
        intercalation_happened = true;
        
        % A. Geometry: Move v2_orphan on top of v1_master
        coordinates(v2_orphan, :) = coordinates(v1_master, :);
        
        % B. Connectivity (Re-wiring):
        for j = 1:length(cells_to_rewire)
            cell_idx = cells_to_rewire(j);
            connec_list = connectivity{cell_idx};
            connec_list(connec_list == v2_orphan) = v1_master;
            connectivity{cell_idx} = connec_list;
        end
        
        % C. Connectivity (Fixing Cell_W):
        % Remove BOTH vertices from Cell_W's list.
        connec_list = connectivity{Cell_W};
        connec_list(connec_list == v1_master) = [];
        connec_list(connec_list == v2_orphan) = [];
        connectivity{Cell_W} = connec_list;
        
        % D. Update Vert-To-Cell Map:
        verttocell{v1_master} = unique([verttocell{v1_master}, cells_to_rewire]);
        vtc_v1_list = verttocell{v1_master};
        vtc_v1_list(vtc_v1_list == Cell_W) = [];
        verttocell{v1_master} = vtc_v1_list;
        verttocell{v2_orphan} = [];
        
        % E. Mark Edge for Removal
        edges_to_remove_idx = [edges_to_remove_idx, i]; %#ok<AGROW>
        
        % F. Set Relaxation Flag
        T1flagVec(v1_master) = 1;
        
        fprintf(1,'--- 3-Cell Intercalation Performed ---\n');
        fprintf(1,'Cell %d fully detached from margin.\n', Cell_W);
        fprintf(1,'Vertex %d merged into %d.\n', v2_orphan, v1_master);
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