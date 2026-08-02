function [P, EdgeData] = getCellPerimeters(nCells,coordinates,connectivity,param, imAugmented)
%% getCellPerimeters: get cell perimeter

P        = zeros(nCells,1);
EdgeData = cell(nCells,1);

for cellID = 1:nCells
	[P(cellID), EdgeData{cellID}] = getPolygonalCellPerimeter(coordinates,connectivity{cellID},param, imAugmented);
end