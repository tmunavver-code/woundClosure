function A = getCellAreas(celldata, param)
%% getCellAreas: provide cell areas (zero for cells with fewer than 3 vertices)

A = zeros(celldata.nCells,1);

for cellID = 1:celldata.nCells
    
    if length(celldata.connec{cellID}) < 3
        A(cellID) = 0; % Invalid remnant cell has zero area
    else
	    A(cellID) = getPolygonalCellArea(celldata,cellID,param);
    end
end
end