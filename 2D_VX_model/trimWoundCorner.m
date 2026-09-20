function connectivity = trimWoundCorner(cellA, connectivity, cornerVertex)
% TRIMWOUNDCORNER Removes a corner vertex touching two wound boundary edges
% from cellA, merging the two boundary edges into one.

vertices = connectivity{cellA};
idx = find(vertices == cornerVertex, 1);
if ~isempty(idx)
    vertices(idx) = [];
    connectivity{cellA} = vertices;
end

end
