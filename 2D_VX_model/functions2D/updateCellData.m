function celldata = updateCellData(x, x_prev, celldata)
%% updateCellData: updates the cell data from current positions of the nodes

celldata.r = reshape(x,2,length(celldata.r))';
celldata.r0 = reshape(x_prev,2,length(celldata.r0))';

celldata.A = getCellAreas(celldata);
celldata.P = getCellPerimeters(celldata);

