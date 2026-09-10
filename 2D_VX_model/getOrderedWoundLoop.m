function [loopVertices, isValid] = getOrderedWoundLoop(woundEdges)
%GETORDEREDWOUNDLOOP Reconstruct a single ordered vertex loop from an
%unordered list of wound-margin edges (celldata.woundEdges, Nx2).
%
% findWoundEdgeCells.m returns the wound margin as an unordered edge
% list. Area/area-derivative code (getPolygonalCellArea.m,
% getAreaDerivative.m) needs an ordered vertex loop, the same way a
% normal cell's celldata.connec{cellID} is ordered.
%
% Assumes the wound margin is a single simple closed loop (every vertex
% has exactly degree 2 in this edge list) -- which is what the "no
% spurious second wound" fix earlier in this project is supposed to
% guarantee. If that assumption doesn't hold (degree != 2 anywhere, or
% the walk doesn't close), isValid=false is returned and the caller
% should skip applying any force this step rather than guess.

isValid = false;
loopVertices = [];

if isempty(woundEdges)
    return;
end

% Build adjacency: vertex -> list of neighbor vertices via wound edges
adj = containers.Map('KeyType','double','ValueType','any');
for i = 1:size(woundEdges,1)
    v1 = woundEdges(i,1);
    v2 = woundEdges(i,2);
    if isKey(adj, v1)
        adj(v1) = [adj(v1), v2];
    else
        adj(v1) = v2;
    end
    if isKey(adj, v2)
        adj(v2) = [adj(v2), v1];
    else
        adj(v2) = v1;
    end
end

allVerts = keys(adj);
for i = 1:length(allVerts)
    if numel(adj(allVerts{i})) ~= 2
        % Not a simple closed loop -- bail out safely.
        return;
    end
end

nVerts = length(allVerts);
startVert = woundEdges(1,1);
loopVertices = zeros(1, nVerts);
loopVertices(1) = startVert;
neighbors = adj(startVert);
prevVert = startVert;
currVert = neighbors(1);

for k = 2:nVerts
    loopVertices(k) = currVert;
    nbrs = adj(currVert);
    if nbrs(1) == prevVert
        nextVert = nbrs(2);
    else
        nextVert = nbrs(1);
    end
    prevVert = currVert;
    currVert = nextVert;
end

% The walk should return exactly to the start.
if currVert == startVert
    isValid = true;
else
    isValid = false;
    loopVertices = [];
end

end
