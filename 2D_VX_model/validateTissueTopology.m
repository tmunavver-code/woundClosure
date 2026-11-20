function [isValid, reason, maxVertexDegree] = validateTissueTopology(celldata)
% VALIDATETISSUETOPOLOGY Checks connectivity and vertex-to-cell consistency.

isValid = false;
nVertices = size(celldata.r, 1);

if any(~isfinite(celldata.r), 'all')
    reason = 'nonfinite_coordinates';
    maxVertexDegree = NaN;
    return;
end

for cellID = 1:celldata.nCells
    vertices = celldata.connec{cellID};
    if isempty(vertices)
        continue;
    end
    if numel(vertices) < 3 || numel(unique(vertices)) < 3
        reason = 'degenerate_cell';
        maxVertexDegree = NaN;
        return;
    end
    if any(vertices < 1) || any(vertices > nVertices) || ...
            any(vertices ~= round(vertices))
        reason = 'invalid_vertex_index';
        maxVertexDegree = NaN;
        return;
    end
    if numel(unique(vertices)) ~= numel(vertices)
        reason = 'repeated_vertex_in_cell';
        maxVertexDegree = NaN;
        return;
    end
end

expectedMap = getVertexToCellIDMap(nVertices, celldata.nCells, celldata.connec);
for vertexID = 1:nVertices
    expected = sort(unique(expectedMap{vertexID}));
    actual = sort(unique(celldata.verttocell{vertexID}));
    if ~isequal(expected, actual)
        reason = 'vertex_to_cell_mismatch';
        maxVertexDegree = NaN;
        return;
    end
end

maxVertexDegree = max(cellfun(@(x) numel(unique(x)), expectedMap));
isValid = true;
reason = 'ok';

end
