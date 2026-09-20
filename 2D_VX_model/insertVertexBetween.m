function [vertsOut, ok] = insertVertexBetween(vertsIn, vA, vB, newV)
% INSERTVERTEXBETWEEN Inserts newV strictly between vA and vB in a cyclic
% vertex list, returning ok=false if vA and vB are not cyclically adjacent.

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
