function [celldata, param] = refreshWoundMargin(celldata, param)
%REFRESHWOUNDMARGIN Recompute the wound margin (free-edge set and the list
%of cells bordering it) from the CURRENT connectivity.
%
% MUST be called after ANY topology change -- bulk T1
% (checkT1transitions.m), cell division (checkCellDivision.m), and
% wound-margin turnover (checkWoundIntercalations.m) all rewire
% connectivity and can change which cells border the wound.
%
% Before this existed, celldata.woundEdges / param.cellIDtoContract were
% only refreshed inside checkWoundIntercalations.m on a successful
% intercalation, so they went stale after every bulk T1 and every
% division. Measured consequences of that staleness in a 1000-step run:
%   - 477 of 487 wound-intercalation candidates that passed the
%     force-based Bell's-law gate were then discarded because the cell
%     actually bordering the wound edge was not in the stale
%     cellIDtoContract list ("cellA ambiguous"), so margin turnover
%     essentially stopped and the margin cell count froze.
%   - The stale woundEdges set stopped matching the live connectivity,
%     which is what made getOrderedWoundLoop fail and the wound-area
%     trace go NaN partway through a run.
%   - param.cellIDtoContract also selects which cells get the
%     wound-margin contractility boost and purse-string tension, so a
%     stale list means those active forces were being applied to the
%     wrong cells -- a physics error, not just bookkeeping.

[woundEdgeCells, woundEdges] = findWoundEdgeCells(celldata);

celldata.woundEdges = woundEdges;
celldata.cellIDtoContract = woundEdgeCells;
param.cellIDtoContract = woundEdgeCells;

end
