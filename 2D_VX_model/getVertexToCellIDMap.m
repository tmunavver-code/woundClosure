function vertocellmap = getVertexToCellIDMap(Nvertices, nCells,connectivity)
% provides vertex to cell number map for later use

vertocellmap = cell(Nvertices,1);

for cellID = 1:nCells
	verticesmat = connectivity{cellID};
	for vertmatind = 1:length(verticesmat)
		vertID = verticesmat(vertmatind);
		vertocellmap{vertID} = [vertocellmap{vertID}, cellID];
	end
end