function dPeri = getPerimeterDerivative(vertexcoords,currVert,nextVert,prevVert)
% get perimeter derivative

RcurrVert = vertexcoords(currVert,:);
RprevVert = vertexcoords(prevVert,:);
RnextVert = vertexcoords(nextVert,:);

dPeri = (RcurrVert - RnextVert)/norm((RcurrVert - RnextVert)) + ...
	(RcurrVert - RprevVert)/norm((RcurrVert - RprevVert));