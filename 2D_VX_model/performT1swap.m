function [coordinates, connectivity, edgedata,verttocell] = performT1swap(cellID1, cellID2, cellID3, cellID4, coordinates,connectivity,edgedata,verttocell, node1, node2, param, imAugmented, T1_TOL_triggered)
%
% MODIFIED:
% 1. Function now accepts 'T1_TOL_triggered' as the last argument.
% 2. Replaced 'param.T1_TOL' with 'T1_TOL_triggered' for vertex bounce.
%

Pos1 = coordinates(node1,:);
Pos2 = coordinates(node2,:);

% First we separate the two nodes further away from the threshold
vec = Pos2 - Pos1;
newPos1 = Pos1 - 0.52 * T1_TOL_triggered * vec/norm(vec); % <-- MODIFIED
newPos2 = Pos2 + 0.52 * T1_TOL_triggered * vec/norm(vec); % <-- MODIFIED
coordinates(node1,:) = newPos1;
coordinates(node2,:) = newPos2;
normNewvec = norm(newPos2 - newPos1);

% Cell 1 (remove node2)
if cellID1 ~= -1
    vertices = connectivity{cellID1};
    node2ind = find(vertices == node2);
    vertices(node2ind) = [];
    connectivity{cellID1} = vertices;
    
    % Removing the edge between node1 and node2 and recalculating edge
    % distance
    node1ind = find(vertices == node1);
    edgedata{cellID1}(node1ind) = [];
    [~,edgedata{cellID1}] = getPolygonalCellPerimeter(coordinates,connectivity{cellID1},param, imAugmented);
    
    % Remove cell1 from node 2
    vtc = verttocell{node2};
    cellID1ind = find(vtc == cellID1);
    vtc(cellID1ind) = [];
    verttocell{node2} = vtc;
end

% Cell 2 (remove node1)
if cellID2 ~= -1
    vertices = connectivity{cellID2};
    node1ind = find(vertices == node1);
    vertices(node1ind) = [];
    connectivity{cellID2} = vertices;
    
    % Removing the edge between node1 and node2 and recalculating edge
    % distance
    node2ind = find(vertices == node2);
    edgedata{cellID2}(node2ind) = [];
    [~,edgedata{cellID2}] = getPolygonalCellPerimeter(coordinates,connectivity{cellID2},param, imAugmented);

    % Remove cell2 from node 1
    vtc = verttocell{node1};
    cellID2ind = find(vtc == cellID2);
    vtc(cellID2ind) = [];
    verttocell{node1} = vtc;
end

% Cell 3 (add node2 before node1)
if cellID3 ~= -1
    vertices = connectivity{cellID3};
    node1ind = find(vertices == node1);
    vertices = insertelem(node2, node1ind , vertices);
    connectivity{cellID3} = vertices;
    
    % Adding the new edge and recalculating edge distance
    EdgeLength = norm(Pos2 - Pos1);
    EdgeVec = edgedata{cellID3};
    EdgeVec = insertelem(EdgeLength,node1ind,EdgeVec);
    edgedata{cellID3} = EdgeVec;
    [~,edgedata{cellID3}] = getPolygonalCellPerimeter(coordinates,connectivity{cellID3},param, imAugmented);
   
    % add cell3 to node 2
    verttocell{node2} = [verttocell{node2}, cellID3];
    
    
end

% Cell 4 (add node1 before node2)
if cellID4 ~= -1
    vertices = connectivity{cellID4};
    node2ind = find(vertices == node2);
    vertices = insertelem(node1, node2ind , vertices);
    connectivity{cellID4} = vertices;
    
    % Adding the new edge and recalculating edge distance
    EdgeLength = norm(Pos2 - Pos1);
    EdgeVec = edgedata{cellID4};
    EdgeVec = insertelem(EdgeLength,node2ind,EdgeVec);
    edgedata{cellID4} = EdgeVec;
    [~,edgedata{cellID4}] = getPolygonalCellPerimeter(coordinates,connectivity{cellID4},param, imAugmented);
    
    % add cell4 to node 1
    verttocell{node1} = [verttocell{node1}, cellID4];
    

end
end