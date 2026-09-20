function [celldata, success] = performCellDivision(cellID, celldata, param)
% PERFORMCELLDIVISION Splits cellID into two daughter cells along its short axis
% passing through the centroid. Returns success=false if division cannot be formed cleanly.

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

% Signed distance of each vertex along long axis relative to centroid
signedDist = zeros(n,1);
for k = 1:n
    signedDist(k) = dot(vertexcoords(vertices(k),:) - centroid, longAxisDir);
end

crossings = []; % edge indices where sign changes
for k = 1:n
    k2 = mod(k, n) + 1;
    if signedDist(k) == 0 || signedDist(k2) == 0
        return; % degenerate (vertex exactly on line)
    end
    if sign(signedDist(k)) ~= sign(signedDist(k2))
        crossings(end+1) = k; %#ok<AGROW>
    end
end

if length(crossings) ~= 2
    return; % non-convex or degenerate
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

% Wrap new vertex coordinates back into periodic box
p1 = mod(p1, [Lx, Ly]);
p2 = mod(p2, [Lx, Ly]);

% Prospective IDs for the two new vertices
newV1 = celldata.nMasterVertices + 1;
newV2 = celldata.nMasterVertices + 2;

% Validate neighbor rewiring before mutating connectivity
neighborsAtEdge1 = setdiff(intersect(celldata.verttocell{v_i}, celldata.verttocell{v_iNext}), cellID);
neighborsAtEdge2 = setdiff(intersect(celldata.verttocell{v_j}, celldata.verttocell{v_jNext}), cellID);

cellN1 = []; connecN1 = [];
if ~isempty(neighborsAtEdge1)
    cellN1 = neighborsAtEdge1(1);
    [connecN1, okN1] = insertVertexBetween(celldata.connec{cellN1}, v_i, v_iNext, newV1);
    if ~okN1
        return;
    end
end

cellN2 = []; connecN2 = [];
if ~isempty(neighborsAtEdge2)
    cellN2 = neighborsAtEdge2(1);
    [connecN2, okN2] = insertVertexBetween(celldata.connec{cellN2}, v_j, v_jNext, newV2);
    if ~okN2
        return;
    end
end

% Build daughter vertex lists by walking original loop
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
    return;
end

% Commit new vertices and updated connectivity
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

% Update vertex-to-cell maps for daughter2
for k = 2:length(daughter2)-1
    v = daughter2(k);
    vtc = celldata.verttocell{v};
    vtc(vtc == cellID) = newCellID;
    celldata.verttocell{v} = vtc;
end

% Commit daughter cells
celldata.connec{cellID} = daughter1;
celldata.connec{newCellID} = daughter2;
celldata.nCells = newCellID;

% Update vertex-to-cell maps for new vertices
vtcNew1 = [cellID, newCellID];
if ~isempty(cellN1), vtcNew1 = [vtcNew1, cellN1]; end
celldata.verttocell{newV1} = vtcNew1;

vtcNew2 = [cellID, newCellID];
if ~isempty(cellN2), vtcNew2 = [vtcNew2, cellN2]; end
celldata.verttocell{newV2} = vtcNew2;

% Initialize daughter cell properties
celldata.A(newCellID, 1) = 0;
celldata.Ainit(newCellID, 1) = celldata.Ainit(cellID);
celldata.P(newCellID, 1) = 0;

success = true;

end
