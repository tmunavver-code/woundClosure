function [P, EdgeData] = getCellPerimeters(nCells,coordinates,connectivity,param, imAugmented)
%% getCellPerimeters: get cell perimeters (zero for cells with fewer than 3 vertices)

P        = zeros(nCells,1);
EdgeData = cell(nCells,1);

for cellID = 1:nCells
    
    vertices = connectivity{cellID};
    
    if length(vertices) < 3
        P(cellID) = 0;
        EdgeData{cellID} = [];
    else
	    [P(cellID), EdgeData{cellID}] = getPolygonalCellPerimeter(coordinates,vertices,param, imAugmented);
    end
end
end