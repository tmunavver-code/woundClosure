function verttocell = getVertexToCellIDMap(nVertices, nCells, connec)
    % GETVERTEXTOCELLIDMAP Creates a mapping from vertex ID to the cells that share it.
    verttocell = cell(nVertices, 1);
    for cellID = 1:nCells
        vertices = connec{cellID};
        for i = 1:length(vertices)
            v = vertices(i);
            if v > 0 && v <= nVertices
                verttocell{v}(end+1) = cellID;
            end
        end
    end
    for i = 1:nVertices
        verttocell{i} = unique(verttocell{i});
    end
end
