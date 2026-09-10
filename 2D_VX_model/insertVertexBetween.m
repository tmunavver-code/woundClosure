function [vertsOut, ok] = insertVertexBetween(vertsIn, vA, vB, newV)
%INSERTVERTEXBETWEEN Insert newV strictly between vA and vB in a cyclic
%vertex list, whichever order they appear in.
%
% Returns ok=false (and vertsIn unchanged) if vA and vB are not both
% present and cyclically adjacent -- the caller must then abort whatever
% topology change it was attempting rather than commit a corrupted mesh.
%
% Needed because a plain "insert before index of vA" is order-dependent:
% adjacent cells traverse a shared edge in opposite winding, so inserting
% before vA lands the new vertex OUTSIDE the vA-vB edge whenever the
% neighbour happens to list them in the same order as the dividing cell.
% That leaves the neighbour's edge unmatched against the daughter cell's
% two half-edges, and findWoundEdgeCells.m then reports both as free
% edges -- i.e. a spurious wound boundary appears in the bulk.

vertsOut = vertsIn;
ok = false;

n = numel(vertsIn);
iA = find(vertsIn == vA, 1);
iB = find(vertsIn == vB, 1);
if isempty(iA) || isempty(iB)
    return;
end

if mod(iA, n) + 1 == iB
    % vB immediately follows vA -- insert between them, i.e. before vB.
    vertsOut = insertelem(newV, iB, vertsIn);
    ok = true;
elseif mod(iB, n) + 1 == iA
    % vA immediately follows vB -- insert before vA.
    vertsOut = insertelem(newV, iA, vertsIn);
    ok = true;
else
    % Not cyclically adjacent -- corrupted or unexpected topology.
    return;
end

end
