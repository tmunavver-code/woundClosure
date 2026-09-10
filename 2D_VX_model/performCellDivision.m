function [celldata, success] = performCellDivision(cellID, celldata, param)
%PERFORMCELLDIVISION Split cellID into two daughter cells along its short
%axis (cleavage line through the centroid, perpendicular to the cell's
%long axis) -- matching Tetley et al. 2019's stated division geometry
%("division of elongated cells along their short axis"), but triggered by
%a force-based area criterion (see checkCellDivision.m) rather than their
%fixed-rate cell-cycle clock.
%
% Returns success=false (celldata unchanged) if the short-axis line does
% not cross exactly two edges of the polygon (e.g. a very non-convex
% cell) -- skip rather than guess, consistent with the rest of this
% codebase's topology-change guards.

success = false;
Lx = param.Lx;
Ly = param.Ly;

vertices = celldata.connec{cellID};
n = length(vertices);
if n < 4
    return; % need at least 4 vertices to split into two valid (>=3-vertex) cells
end

vertexcoords = modifyVerticesForPeriodicity(celldata.r, vertices, Lx, Ly);
[longAxisDir, ~, centroid] = getCellPrincipalAxis(vertexcoords, vertices);

% Signed distance of each vertex along the long axis, relative to the
% centroid -- the cleavage line is where this changes sign (it is the
% line through the centroid perpendicular to the long axis, i.e. along
% the short axis).
signedDist = zeros(n,1);
for k = 1:n
    signedDist(k) = dot(vertexcoords(vertices(k),:) - centroid, longAxisDir);
end

crossings = []; % edge indices (k means edge vertices(k)-vertices(k+1)) where sign changes
for k = 1:n
    k2 = mod(k, n) + 1;
    if signedDist(k) == 0 || signedDist(k2) == 0
        return; % degenerate (vertex exactly on the line) -- skip rather than special-case
    end
    if sign(signedDist(k)) ~= sign(signedDist(k2))
        crossings(end+1) = k; %#ok<AGROW>
    end
end

if length(crossings) ~= 2
    return; % non-convex or degenerate -- skip
end

k1 = crossings(1);
k2 = crossings(2);
k1n = mod(k1, n) + 1;
k2n = mod(k2, n) + 1;

v_i  = vertices(k1);  v_iNext = vertices(k1n);
v_j  = vertices(k2);  v_jNext = vertices(k2n);

s1 = signedDist(k1) / (signedDist(k1) - signedDist(k1n));
p1 = vertexcoords(v_i,:) + s1 * (vertexcoords(v_iNext,:) - vertexcoords(v_i,:));
s2 = signedDist(k2) / (signedDist(k2) - signedDist(k2n));
p2 = vertexcoords(v_j,:) + s2 * (vertexcoords(v_jNext,:) - vertexcoords(v_j,:));

% Wrap the new points back into the periodic box (vertexcoords above was
% unwrapped relative to this cell's own reference vertex).
p1 = mod(p1, [Lx, Ly]);
p2 = mod(p2, [Lx, Ly]);

% Prospective IDs for the two new vertices (not committed yet -- the
% neighbour rewiring below is validated first, and the whole division is
% abandoned if it cannot be done cleanly).
newV1 = celldata.nMasterVertices + 1;
newV2 = celldata.nMasterVertices + 2;

% Validate the neighbour rewiring BEFORE mutating anything. Each of the
% two split edges may be shared with a neighbouring cell, which must gain
% the corresponding new vertex at exactly the right place in its own
% (opposite-winding) vertex list -- otherwise its edge no longer matches
% the daughter cells' half-edges and findWoundEdgeCells.m reports
% spurious free edges. An edge with no neighbour (i.e. on the wound
% boundary) is fine: it simply becomes two free edges.
neighborsAtEdge1 = setdiff(intersect(celldata.verttocell{v_i}, celldata.verttocell{v_iNext}), cellID);
neighborsAtEdge2 = setdiff(intersect(celldata.verttocell{v_j}, celldata.verttocell{v_jNext}), cellID);

cellN1 = []; connecN1 = [];
if ~isempty(neighborsAtEdge1)
    cellN1 = neighborsAtEdge1(1);
    [connecN1, okN1] = insertVertexBetween(celldata.connec{cellN1}, v_i, v_iNext, newV1);
    if ~okN1
        return; % cannot rewire cleanly -- abandon this division
    end
end

cellN2 = []; connecN2 = [];
if ~isempty(neighborsAtEdge2)
    cellN2 = neighborsAtEdge2(1);
    [connecN2, okN2] = insertVertexBetween(celldata.connec{cellN2}, v_j, v_jNext, newV2);
    if ~okN2
        return; % cannot rewire cleanly -- abandon this division
    end
end

% Build the two daughter vertex lists by walking the original loop,
% inserting newV1 between v_i/v_iNext and newV2 between v_j/v_jNext.
daughter1 = newV1;
idx = k1n;
while idx ~= k2n
    daughter1(end+1) = vertices(idx); %#ok<AGROW>
    idx = mod(idx, n) + 1;
end
daughter1(end+1) = newV2;

daughter2 = newV2;
idx = k2n;
while idx ~= k1n
    daughter2(end+1) = vertices(idx); %#ok<AGROW>
    idx = mod(idx, n) + 1;
end
daughter2(end+1) = newV1;

if length(daughter1) < 3 || length(daughter2) < 3
    return; % would create a degenerate daughter -- skip
end

% All validation passed -- only now commit anything.
celldata.r(newV1, :) = p1;
celldata.r(newV2, :) = p2;
celldata.nMasterVertices = celldata.nMasterVertices + 2;
if ~isempty(cellN1)
    celldata.connec{cellN1} = connecN1;
end
if ~isempty(cellN2)
    celldata.connec{cellN2} = connecN2;
end

newCellID = celldata.nCells + 1;

% Re-point every vertex that ends up in daughter2 (except the two brand
% new ones, handled separately below) from cellID to newCellID.
for k = 2:length(daughter2)-1
    v = daughter2(k);
    vtc = celldata.verttocell{v};
    vtc(vtc == cellID) = newCellID;
    celldata.verttocell{v} = vtc;
end

% (Neighbour rewiring was validated and committed above, before any other
% mutation -- see insertVertexBetween.m for why the placement has to be
% strictly between the two shared vertices.)

% Commit the two daughter cells (reuse cellID for daughter1, append
% daughter2 as a new cell) and the vertex-to-cell map for the two new
% vertices.
celldata.connec{cellID} = daughter1;
celldata.connec{newCellID} = daughter2;
celldata.nCells = newCellID;

% newV1 lies on the old v_i-v_iNext edge, which daughter1, daughter2 and
% (if present) the neighbour cellN1 all now touch; likewise newV2.
vtcNew1 = [cellID, newCellID];
if ~isempty(cellN1), vtcNew1 = [vtcNew1, cellN1]; end
celldata.verttocell{newV1} = vtcNew1;

vtcNew2 = [cellID, newCellID];
if ~isempty(cellN2), vtcNew2 = [vtcNew2, cellN2]; end
celldata.verttocell{newV2} = vtcNew2;

% Daughter cells inherit the parent's area target bookkeeping (A, Ainit
% get recomputed by the caller via getCellAreas immediately after, same
% pattern as every other topology change in this codebase).
celldata.A(newCellID, 1) = 0;
celldata.Ainit(newCellID, 1) = celldata.Ainit(cellID);
celldata.P(newCellID, 1) = 0;

success = true;

end
