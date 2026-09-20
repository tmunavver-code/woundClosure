function [coordinates, connectivity, verttocell] = performWoundMarginT1(cellA, cellB, cellC, coordinates, connectivity, verttocell, v1_master, v2_orphan, param)
% PERFORMWOUNDMARGINT1 Executes a T1-style neighbor swap at the wound margin,
% treating the wound as the 4th cell in the tetrad (cellA loses margin edge,
% cellB and cellC gain vertices).

Pos1 = coordinates(v1_master,:);
Pos2 = coordinates(v2_orphan,:);
vec = Pos2 - Pos1;
Lsep = norm(vec);
if Lsep < 1e-10
    % Nudge coincident vertices apart along arbitrary direction
    vec = [1, 0];
    Lsep = 1;
end
% Push vertices apart using separation scale proportional to T1_TOL_wound
sepScale = 0.52 * param.T1_TOL_wound;
coordinates(v1_master,:) = Pos1 - sepScale * vec / Lsep;
coordinates(v2_orphan,:) = Pos2 + sepScale * vec / Lsep;

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
