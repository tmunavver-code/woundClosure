function [Vsubset, Csubset] = getCellDataWithinDomain(Ncells, V, C, closedPolygonsID)
%% getCellDataWithinDomain: Extract Vertices and Cell data
% V is the list of vertices, and C is the cell connectivity matrix
% closedPolygonsID contains the list of cellIDs that are required in the
% output voronoi diagram

%% Process voronoi input to keep only the relevant vertices
% those that appear in closedPolygonsID
Vsubset = [];
Csubset = cell(Ncells, 1);
Nvertsubset = 0;
addvertflag = 0;

for cellID = 1:Ncells
    cellIDglobal = closedPolygonsID(cellID);
    vertmat      = C{cellIDglobal};
    Nvertices    = length(vertmat);
    newvertmat   = [];
    
    for vertind = 1:Nvertices
        vertIDglobal = vertmat(vertind);
        Vx           = V(vertIDglobal,1);
        Vy           = V(vertIDglobal,2);
        
        % check if this vertex already exists in Vsubset
        if Nvertsubset > 0
            diff_vert = sqrt((Vsubset(:,1) - Vx).^2 + (Vsubset(:,2) - Vy).^2) < 1E-06;
            vertID_match = find(diff_vert);
            if ~isempty(vertID_match)
                vertID = vertID_match;
            else
                addvertflag = 1;
            end
        else
            addvertflag = 1; % For the first entry
        end
        
        % Add new entry if necessary
        if addvertflag == 1
            Nvertsubset = Nvertsubset + 1;
            vertID      = Nvertsubset;
            Vsubset     = [Vsubset; Vx, Vy];
            addvertflag = 0;
        end
        
        % update newvertmat
        newvertmat = [newvertmat, vertID];
    end
    
    % Update cell entry
    Csubset{cellID} = newvertmat;
end