function [Pcell, EdgesCell] = getPolygonalCellPerimeter(coordinates,connectivityVec,param, imAugmented)

vertices          = connectivityVec;
vertexcoordsold   = coordinates;
Lx                = param.Lx;
Ly                = param.Ly;
Pcell             = 0;
EdgesCell         = zeros(size(vertices));

% Modify verecoords acocuntin for peridicity
if imAugmented == 0
    vertexcoords = modifyVerticesForPeriodicity(vertexcoordsold,vertices,Lx,Ly);
else
    vertexcoords = vertexcoordsold;
end

for j = 1:length(vertices)
	if j < length(vertices)
		p1 = vertices(j);
		p2 = vertices(j+1);
	else % Final edge
		p1 = vertices(j);
		p2 = vertices(1);
	end
	
	vecp1 = [vertexcoords(p1,1),vertexcoords(p1,2),0.0];
	vecp2 = [vertexcoords(p2,1),vertexcoords(p2,2),0.0];
	
	EdgesCell(j) = norm(vecp2 - vecp1);
	Pcell        = Pcell + EdgesCell(j);
end