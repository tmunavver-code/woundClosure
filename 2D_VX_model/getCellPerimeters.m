function [P, EdgeData] = getCellPerimeters(nCells,coordinates,connectivity,param, imAugmented)
%% getCellPerimeters: get cell perimeter
%
% MODIFIED (Stable V2):
% - Added a guard clause to skip perimeter calculations for invalid cells
%   (with < 3 vertices), setting their perimeter to 0.
%

P        = zeros(nCells,1);
EdgeData = cell(nCells,1);

for cellID = 1:nCells
    
    vertices = connectivity{cellID};
    
    % --- NEW GUARD CLAUSE ---
    if length(vertices) < 3
        P(cellID) = 0;
        EdgeData{cellID} = [];
    else
    % --- END NEW ---
	    [P(cellID), EdgeData{cellID}] = getPolygonalCellPerimeter(coordinates,vertices,param, imAugmented);
    end
end
end