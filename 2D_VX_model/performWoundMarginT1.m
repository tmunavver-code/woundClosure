function [coordinates, connectivity, verttocell] = performWoundMarginT1(cellA, cellB, cellC, coordinates, connectivity, verttocell, v1_master, v2_orphan, param)
% PERFORMWOUNDMARGINT1 Executes a T1-style neighbor swap at the wound margin,
% treating the wound as the 4th cell in the tetrad (cellA loses margin edge,
% cellB and cellC gain vertices).

Pos1 = coordinates(v1_master,:);
Pos2 = coordinates(v2_orphan,:);
oldEdge = Pos2 - Pos1;
oldEdge(1) = oldEdge(1) - param.Lx * round(oldEdge(1) / param.Lx);
oldEdge(2) = oldEdge(2) - param.Ly * round(oldEdge(2) / param.Ly);
oldLength = norm(oldEdge);
if oldLength < 1e-10
    oldDirection = [1, 0];
else
    oldDirection = oldEdge / oldLength;
end

% A continuous T1 passes through a zero-length fourfold vertex.  Place the
% replacement edge at a small length perpendicular to the disappearing
% edge, rather than expanding the old edge by a finite distance.
midpoint = Pos1 + 0.5 * oldEdge;
newDirection = [-oldDirection(2), oldDirection(1)];
if isfield(param, 'woundT1NewEdgeLength')
    newLength = param.woundT1NewEdgeLength;
else
    newLength = 0.01 * param.T1_TOL_wound;
end
coordinates(v1_master,:) = midpoint - 0.5 * newLength * newDirection;
coordinates(v2_orphan,:) = midpoint + 0.5 * newLength * newDirection;
coordinates(:,1) = mod(coordinates(:,1), param.Lx);
coordinates(:,2) = mod(coordinates(:,2), param.Ly);

% Cell A: remove v2_orphan (keeps v1_master) -- pushed off the margin.
verticesA = connectivity{cellA};
idxA = find(verticesA == v2_orphan, 1);
if ~isempty(idxA)
    verticesA(idxA) = [];
    connectivity{cellA} = verticesA;
end
vtc = verttocell{v2_orphan};
vtc(vtc == cellA) = [];
verttocell{v2_orphan} = vtc;

% Cell B: insert v2_orphan immediately before v1_master (gains a vertex).
verticesB = connectivity{cellB};
idxB = find(verticesB == v1_master, 1);
if ~isempty(idxB)
    verticesB = insertelem(v2_orphan, idxB, verticesB);
    connectivity{cellB} = verticesB;
end
verttocell{v2_orphan} = [verttocell{v2_orphan}, cellB];

% Cell C: insert v1_master immediately before v2_orphan (gains a vertex).
verticesC = connectivity{cellC};
idxC = find(verticesC == v2_orphan, 1);
if ~isempty(idxC)
    verticesC = insertelem(v1_master, idxC, verticesC);
    connectivity{cellC} = verticesC;
end
verttocell{v1_master} = [verttocell{v1_master}, cellC];

end
