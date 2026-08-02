function [celldata, T1flagVec, intercalation_happened] = checkWoundIntercalations(celldata, param, T1flagVec)
%
% Implements a 3-cell intercalation (topologically a T2 extrusion)
% at the wound margin using a "Vertex Merge & Orphan" method.
%
% THIS IS VERSION 1: It has the bug where Cell_W remains
% attached to the master vertex (v1_master).
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
    % Calculate edge length with Minimum Image Convention (MIC)
    vec = coordinates(v2_orphan,:) - coordinates(v1_master,:);
    vec(1) = vec(1) - param.Lx * round(vec(1) / param.Lx);
    vec(2) = vec(2) - param.Ly * round(vec(2) / param.Ly);
    L = norm(vec);
    
    if L < param.L_intercalation_thresh
        
        % --- 2. IDENTIFY ALL PLAYERS ---
        
        % Get all cells touching the two vertices
        cells_at_v1 = verttocell{v1_master};
        cells_at_v2 = verttocell{v2_orphan};
        
        % Find the cell that touches BOTH (this is the one to extrude)
        Cell_W_candidates = intersect(cells_at_v1, cells_at_v2);
        
        % Of these, find the one that is *also* a wound cell
        Cell_W = intersect(Cell_W_candidates, woundCells);
        
        % GUARD CLAUSE: This logic is *only* for a 3-cell junction.
        if length(Cell_W) ~= 1
            continue;
        end
        
        % Find cells to be RE-WIRED. These are cells that touch v2_orphan
        % but NOT v1_master.
        cells_to_rewire = setdiff(cells_at_v2, cells_at_v1);
        
        % --- 3. PERFORM VERTEX MERGE & ORPHAN ---
        intercalation_happened = true;
        
        % A. Geometry: Move v2_orphan on top of v1_master
        coordinates(v2_orphan, :) = coordinates(v1_master, :);
        
        % B. Connectivity (Re-wiring):
        % For every cell that touched *only* v2, re-wire it to v1.
        for j = 1:length(cells_to_rewire)
            cell_idx = cells_to_rewire(j);
            connec_list = connectivity{cell_idx};
            connec_list(connec_list == v2_orphan) = v1_master;
            connectivity{cell_idx} = connec_list;
        end
        
        % C. Connectivity (Fixing Cell_W):
        % Remove the orphaned vertex v2_orphan from Cell_W's list.
        % *** THIS IS THE BUGGY PART ***
        % It *should* also remove v1_master but doesn't.
        connec_list = connectivity{Cell_W};
        connec_list(connec_list == v2_orphan) = [];
        connectivity{Cell_W} = connec_list;
        
        % D. Update Vert-To-Cell Map:
        % Move all re-wired cells from v2_orphan to v1_master
        verttocell{v1_master} = unique([verttocell{v1_master}, cells_to_rewire]);
        % Orphan v2_orphan (it no longer belongs to any cell)
        verttocell{v2_orphan} = [];
        
        % E. Mark Edge for Removal
        edges_to_remove_idx = [edges_to_remove_idx, i]; %#ok<AGROW>
        
        % F. Set Relaxation Flag
        T1flagVec(v1_master) = 1;
        
        fprintf(1,'--- 3-Cell Intercalation Performed ---\n');
        fprintf(1,'Cell %d extruded from margin.\n', Cell_W);
        fprintf(1,'Vertex %d merged into %d.\n', v2_orphan, v1_master);
    end
end

% --- 4. CLEAN UP & UPDATE STATE ---
if intercalation_happened
    % Commit connectivity & coordinate changes
    celldata.r = coordinates;
    celldata.connec = connectivity;
    celldata.verttocell = verttocell;
    
    % The wound margin is now different. We MUST find the new
    % wound edges and wound cells.
    [newWoundCells, newWoundEdges] = findWoundEdgeCells(celldata);
    
    celldata.woundEdges = newWoundEdges;
    param.cellIDtoContract = newWoundCells; % Update list for plotting
    celldata.cellIDtoContract = newWoundCells; % Also store in celldata
end

end