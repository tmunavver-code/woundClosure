function connectivity = trimWoundCorner(cellA, connectivity, cornerVertex)
%TRIMWOUNDCORNER Remove a "corner" vertex from a wound-margin cell.
%
% A corner vertex is one where BOTH of its edges are wound-boundary
% edges -- i.e. cellA is the only real cell touching it (both of its
% other neighboring faces are already the wound). A standard 4-cell T1
% handoff (performWoundMarginT1.m) has no partner cell to give a vertex
% to in this configuration, which otherwise freezes margin turnover
% entirely once the shrinking wound makes corners common (see the
% surrounding checkWoundIntercalations.m comment for the empirical
% symptom this fixes).
%
% Removing the corner vertex from cellA merges its two wound-boundary
% edges into one longer edge -- cellA loses one vertex (the "spike" is
% trimmed) but otherwise stays a completely normal, intact cell, exactly
% like a standard T1: no cell is destroyed by this operation.

vertices = connectivity{cellA};
idx = find(vertices == cornerVertex, 1);
if ~isempty(idx)
    vertices(idx) = [];
    connectivity{cellA} = vertices;
end

end
