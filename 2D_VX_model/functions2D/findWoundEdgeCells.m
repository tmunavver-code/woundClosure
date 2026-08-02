function [woundEdgeCells, woundEdges] = findWoundEdgeCells(celldata)
%
% High-performance replacement for findWoundEdgeCells.
% Uses a containers.Map (hash map) to find free edges in O(N) time
% instead of the original O(N^2) nested loop.
%

% A hash map to store edges.
% Key: string (e.g., '10_25')
% Value: {cellID, [v1, v2]}
edgeMap = containers.Map('KeyType', 'char', 'ValueType', 'any');

for i = 1:length(celldata.connec)
    vertices = celldata.connec{i};
    n = length(vertices);
    
    if n < 2 % Skip cells that have been reduced to a point
        continue;
    end
    
    for j = 1:n
        v1 = vertices(j);
        if j < n
            v2 = vertices(j+1);
        else
            v2 = vertices(1);  % wrap around to first vertex
        end
        
        % Create a unique, direction-independent key for the edge
        if v1 < v2
            edgeKey = [num2str(v1), '_', num2str(v2)];
            edge = [v1, v2];
        else
            edgeKey = [num2str(v2), '_', num2str(v1)];
            edge = [v2, v1];
        end
        
        % --- This is the O(N) logic ---
        if isKey(edgeMap, edgeKey)
            % Key already exists, meaning this is an internal edge
            % shared by two cells. Remove it.
            remove(edgeMap, edgeKey);
        else
            % Key does not exist. This is the first time we've seen
            % this edge, so it's a potential wound edge.
            % Store the cellID and the edge vertices.
            edgeMap(edgeKey) = {i, edge};
        end
    end
end

% --- Post-processing ---
% At this point, edgeMap contains *only* the wound edges.
% We just need to format them into the output variables.
allMapValues = values(edgeMap);

if isempty(allMapValues)
    % No wound edges found
    woundEdgeCells = [];
    woundEdges = zeros(0, 2); % Return correct empty matrix size
    return;
end

% Pre-allocate for speed
numWoundEdges = length(allMapValues);
cellList = zeros(numWoundEdges, 1);
woundEdges = zeros(numWoundEdges, 2);

% Loop through the remaining values (wound edges)
for k = 1:numWoundEdges
    valueData = allMapValues{k}; % valueData is {cellID, [v1, v2]}
    cellList(k) = valueData{1};
    woundEdges(k, :) = valueData{2};
end

% Get the unique list of cell IDs
woundEdgeCells = unique(cellList);

end