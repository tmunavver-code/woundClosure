function [celldata] = boundaryVertices(celldata,param)
%BOUNDARYVERTICES Identifying the boundary vertices through the slaves

coordsSlaves = celldata.slave_r;

idxslaveTopBoundary = find(coordsSlaves(:,2) > param.Ly);
idxslaveBottomBoundary = find(coordsSlaves(:,2) < 0);
idxslaveLeftBoundary = find(coordsSlaves(:,1) < 0);
idxslaveRightBoundary = find(coordsSlaves(:,1) > param.Lx);


celldata.idxTopBoundary = celldata.StoMmap(idxslaveBottomBoundary);
celldata.idxBottomBoundary = celldata.StoMmap(idxslaveTopBoundary);
celldata.idxLeftBoundary = celldata.StoMmap(idxslaveRightBoundary);
celldata.idxRightBoundary = celldata.StoMmap(idxslaveLeftBoundary);

celldata.boundaryNodes = [celldata.idxTopBoundary; celldata.idxBottomBoundary;  celldata.idxLeftBoundary; celldata.idxRightBoundary];

end

