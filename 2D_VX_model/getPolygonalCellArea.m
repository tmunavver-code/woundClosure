function Acell = getPolygonalCellArea(celldata,cellID,param)
%% get polygonal area of cellID

vertices       = celldata.connec{cellID};
vertexcoordsold   = celldata.r;
Lx             = param.Lx;
Ly             = param.Ly;
Acell          = 0;

% Modify vertex coords accounting for peridicity
vertexcoords = modifyVerticesForPeriodicity(vertexcoordsold,vertices,Lx,Ly);

% find mean of vertex coordinates
meanCoord = zeros(1,2);
for k = 1:length(vertices)
	meanCoord = meanCoord + vertexcoords(vertices(k),:);
end
meanCoord = meanCoord/length(vertices);

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
	vecp3 = [meanCoord(1),meanCoord(2),0.0];
	
	% 1,2,3 in counterclockwise order
	Acell = Acell + getTriangleArea(vecp1,vecp2,vecp3);
end