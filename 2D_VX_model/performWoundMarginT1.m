function [coordinates, connectivity, verttocell] = performWoundMarginT1(cellA, cellB, cellC, coordinates, connectivity, verttocell, v1_master, v2_orphan, param)
%PERFORMWOUNDMARGINT1 Genuine T1-style neighbor swap at the wound margin,
%treating the wound itself as the 4th "cell" in the tetrad -- matching
%Tetley et al. 2019's own stated method ("the wound itself was treated as
%a cell, meaning an intercalating tetrad was formed of three wound-edge
%cells and the wound") and confirmed against their Supplementary Video 4:
%a cell pushed off the wound margin remains a fully intact, normal cell,
%it does not shrink or vanish.
%
% Mirrors performT1swap.m exactly for the three real cells involved; the
% wound has no real connec/verttocell entry, so there is no "cellID2"
% bookkeeping to do for it -- its new boundary falls out automatically
% once the caller recomputes findWoundEdgeCells afterward.
%
%   cellA: currently borders both v1_master and v2_orphan (touches the
%          wound along this edge). Loses v2_orphan -- pushed off the
%          margin, stays a complete cell.
%   cellB: borders v1_master only. Gains v2_orphan.
%   cellC: borders v2_orphan only. Gains v1_master.

Pos1 = coordinates(v1_master,:);
Pos2 = coordinates(v2_orphan,:);
vec = Pos2 - Pos1;
Lsep = norm(vec);
if Lsep < 1e-10
    % Coincident points (already merged by an older run's state, or a
    % genuinely degenerate edge) -- nudge apart along an arbitrary
    % direction so the swap is still geometrically well defined.
    vec = [1, 0];
    Lsep = 1;
end
% BUG FIX: use the same separation scale bulk T1 uses (T1_TOL_wound),
% NOT L_intercalation_thresh (0.5) -- that's the *trigger* threshold for
% this event, much larger than bulk T1's, and was fine for the old
% near-zero-cost merge mechanism. Reused as the actual geometric
% push-apart distance for a real T1-style flip, it produced a 3-10x
% larger perturbation than bulk T1 ever does, and energy cost scales
% roughly with the square of that -- empirically this caused every
% attempt to cost dE ~5-9, an order of magnitude above what
% T_eff_T2 was calibrated against, so the energy/Metropolis gate
% rejected essentially everything.
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
