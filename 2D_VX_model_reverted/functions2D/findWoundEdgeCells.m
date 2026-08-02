function [woundEdgeCells, woundEdges] = findWoundEdgeCells(celldata)

% Step 1: Gather all edges and store which cell they came from
allEdges = [];     % will store [v1, v2] for each edge
cellIDs  = [];     % matching cell ID for each edge

for i = 1:length(celldata.connec)
    vertices = celldata.connec{i};
    n = length(vertices);
    
    for j = 1:n
        v1 = vertices(j);
        if j < n
            v2 = vertices(j+1);
        else
            v2 = vertices(1);  % wrap around to first vertex
        end
        
        % Make edge direction-independent: sort vertices
        edge = sort([v1, v2]);
        
        % Add to list
        allEdges = [allEdges; edge];     % each row is an edge
        cellIDs  = [cellIDs; i];         % keep track of source cell
    end
end

% Step 2: Find edges that appear only once
numEdges = size(allEdges, 1);
woundEdgeCells = [];

for i = 1:numEdges
    currentEdge = allEdges(i, :);
    count = 0;
    
    for j = 1:numEdges
        if all(allEdges(j, :) == currentEdge)
            count = count + 1;
        end
    end
    
    % If the edge appears only once, it's a wound edge
    if count == 1
        woundEdgeCells = [woundEdgeCells; cellIDs(i)];
    end
end

woundEdges = [];
for i = 1:size(allEdges, 1)
    edge = allEdges(i, :);
    count = sum(all(allEdges == edge, 2));
    if count == 1
        woundEdges = [woundEdges; edge];  % this edge appears only once
    end
end

% Step 3: Keep only unique cell IDs
woundEdgeCells = unique(woundEdgeCells);

end