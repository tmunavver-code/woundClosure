function [ baryCoords ] = getCellBarycenter( cellID, Connectivity, AllverticesCoords )
%GETCELLBARYCENTER Getting the barycenter position of a given cell

VecOfVertexIDinCell = Connectivity{cellID};
CoordsOfAllVerticesInCell = AllverticesCoords(VecOfVertexIDinCell,:);
baryCoords = mean(CoordsOfAllVerticesInCell,1);

end

