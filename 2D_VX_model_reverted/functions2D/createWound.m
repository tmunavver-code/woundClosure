function [celldata, cellsToRemove] = createWound(celldata, woundShape, params, param)
%CREATEWOUND Removes cells from celldata to create a wound of a specific shape.
%   FINAL CORRECTED VERSION:
%   - Calculates cell barycenters correctly for periodic boundary conditions.
%   - Uses periodic distance to select cells for removal.

% Step 1: Calculate the PERODIC geometric center (barycenter) of every cell
numCells = celldata.nCells;
cellCenters = zeros(numCells, 2);
Lx = param.Lx;
Ly = param.Ly;

for i = 1:numCells
    % --- NEW SECTION: Periodic Barycenter Calculation ---
    vertexIDs = celldata.connec{i};
    vertexCoords = celldata.r(vertexIDs, :);
    
    % Use the first vertex as a reference point
    refPoint = vertexCoords(1, :);
    unwrappedCoords = zeros(size(vertexCoords));
    unwrappedCoords(1, :) = refPoint;
    
    % "Unwrap" other vertices relative to the first one
    for j = 2:length(vertexIDs)
        pt = vertexCoords(j, :);
        dx = pt(1) - refPoint(1);
        dy = pt(2) - refPoint(2);
        
        if dx > Lx/2
            pt(1) = pt(1) - Lx;
        elseif dx < -Lx/2
            pt(1) = pt(1) + Lx;
        end
        
        if dy > Ly/2
            pt(2) = pt(2) - Ly;
        elseif dy < -Ly/2
            pt(2) = pt(2) + Ly;
        end
        unwrappedCoords(j, :) = pt;
    end
    
    % Calculate the mean of the UNWRAPPED coordinates
    barycenter = mean(unwrappedCoords, 1);
    
    % Re-wrap the final barycenter coordinate into the box
    cellCenters(i, :) = mod(barycenter, [Lx, Ly]);
    % --- END NEW SECTION ---
end

% Step 2: Identify which cells to remove based on the chosen shape
cellsToRemove = [];
switch woundShape
    case 'circular'
        center = params.center;
        radius_sq = params.radius^2;
        for i = 1:numCells
            dx = abs(cellCenters(i,1) - center(1));
            dy = abs(cellCenters(i,2) - center(2));
            if dx > Lx/2, dx = Lx - dx; end
            if dy > Ly/2, dy = Ly - dy; end
            dist_sq = dx^2 + dy^2;

            if dist_sq < radius_sq
                cellsToRemove = [cellsToRemove, i];
            end
        end

    case 'elliptical'
        center = params.center;
        a_sq = params.x_axis^2;
        b_sq = params.y_axis^2;
        for i = 1:numCells
            dx = abs(cellCenters(i,1) - center(1));
            dy = abs(cellCenters(i,2) - center(2));
            if dx > Lx/2, dx = Lx - dx; end
            if dy > Ly/2, dy = Ly - dy; end
            val = (dx^2 / a_sq) + (dy^2 / b_sq);
            
            if val < 1
                cellsToRemove = [cellsToRemove, i];
            end
        end
        
    otherwise
        warning('Unknown wound shape specified. No wound created.');
        return;
end

if isempty(cellsToRemove)
    warning('No cells found within the specified wound region. Check parameters.');
    return;
end

fprintf('Creating %s wound: Removing %d cells.\n', woundShape, length(cellsToRemove));

% Step 3: Remove the identified cells and their data from celldata
allCellIDs = 1:numCells;
cellsToKeep = setdiff(allCellIDs, cellsToRemove);

celldata.connec = celldata.connec(cellsToKeep);
celldata.A = celldata.A(cellsToKeep);
celldata.Ainit = celldata.Ainit(cellsToKeep);
celldata.P = celldata.P(cellsToKeep);

celldata.nCells = length(cellsToKeep);

celldata.verttocell = getVertexToCellIDMap(celldata.nMasterVertices, celldata.nCells, celldata.connec);

end