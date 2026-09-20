function [celldata, param] = refreshWoundMargin(celldata, param)
% REFRESHWOUNDMARGIN Recomputes wound margin (free edges and boundary cells)
% from the current connectivity after topological events.

[woundEdgeCells, woundEdges] = findWoundEdgeCells(celldata);

celldata.woundEdges = woundEdges;
celldata.cellIDtoContract = woundEdgeCells;
param.cellIDtoContract = woundEdgeCells;

end
