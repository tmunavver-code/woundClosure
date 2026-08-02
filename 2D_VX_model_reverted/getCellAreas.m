function A = getCellAreas(celldata, param)
%% getCellAreas: provide cell areas

A = zeros(celldata.nCells,1);

for cellID = 1:celldata.nCells
	A(cellID) = getPolygonalCellArea(celldata,cellID,param);
end