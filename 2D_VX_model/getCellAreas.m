function A = getCellAreas(celldata, param)
%% getCellAreas: provide cell areas
%
% MODIFIED (Stable V2):
% - Added a guard clause to skip area calculations for invalid cells
%   (with < 3 vertices), setting their area to 0.
%

A = zeros(celldata.nCells,1);

for cellID = 1:celldata.nCells
    
    % --- NEW GUARD CLAUSE ---
    if length(celldata.connec{cellID}) < 3
        A(cellID) = 0; % Invalid remnant cell has zero area
    else
    % --- END NEW ---
	    A(cellID) = getPolygonalCellArea(celldata,cellID,param);
    end
end
end