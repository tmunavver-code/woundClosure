function [cellID1, cellID2, cellID3, cellID4] = getT1cellIDs(node1, node2, verttocell,connectivity)
% Find the cells corresponding to the edge to be flipped to perform the T1
% transition.

% Initialize outputs — guarantees all are assigned even if topology is invalid
cellID1 = -1; cellID2 = -1; cellID3 = -1; cellID4 = -1;

connectedcellsnode1 = verttocell{node1};
connectedcellsnode2 = verttocell{node2};

% Matrix of the cellIDs surrounding the edge
connectedcells = unique([connectedcellsnode1, connectedcellsnode2]);

if length(connectedcells) ~=4
	fprintf(1,'Less then 4 surrounding cells found. Most likely a boundary vertex.\n');
	% do nothing
	cellID1 = -1;
	cellID2 = -1;
	cellID3 = -1;
	cellID4 = -1;
else
	for cellind = 1:length(connectedcells)
		cellID = connectedcells(cellind);
		vertices = connectivity{cellID};

		node1ind = find(vertices == node1);
		node2ind = find(vertices == node2);
		
		node1nextind = node1ind + 1;
		node1prevind = node1ind - 1;
		if node1ind == length(vertices)
			node1nextind = 1;
		elseif node1ind == 1
			node1prevind = length(vertices);
		end
		
		if ~isempty(node1ind)
			if isempty(node2ind)
				cellID3 = cellID;
			else
				if node2ind == node1nextind
					cellID1 = cellID;
				elseif node2ind == node1prevind
					cellID2 = cellID;
				else
					% Corrupted topology: both nodes present but not adjacent
					cellID1 = -1; cellID2 = -1; cellID3 = -1; cellID4 = -1;
					return;
				end
			end
		elseif ~isempty(node2ind)
			cellID4 = cellID;
		else
			% Corrupted topology: cell has neither node
			cellID1 = -1; cellID2 = -1; cellID3 = -1; cellID4 = -1;
			return;
		end
	end
end

