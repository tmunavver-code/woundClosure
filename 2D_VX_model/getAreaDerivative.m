function dA2d    = getAreaDerivative(vertexcoords,vertices,currVert,nextVert,prevVert)
% derivative of cell area with respect to verex position

N_v = length(vertices);

% find mean of vertex coordinates (pericenter coordinate)
Rp = zeros(1,3);
for k = 1:N_v
	Rp(1:2) = Rp(1:2) + vertexcoords(vertices(k),:);
end
Rp = Rp/N_v;

RcurrVert = [vertexcoords(currVert,:),0];
RnextVert = [vertexcoords(nextVert,:),0];
RprevVert = [vertexcoords(prevVert,:),0];

% Current and next vertex
% Area of the triangle
S = 0.5 * cross3d((Rp - RnextVert), (RcurrVert - RnextVert));
modS = norm(S);

dA = 0.5 * cross3d(S/modS, (Rp - RnextVert));% + 0.5/N_v * cross((RcurrVert - RnextVert),S/modS);

% Current and prev vertex
S = 0.5 * cross3d((RcurrVert - RprevVert), (Rp - RprevVert));
modS = norm(S);

dA = dA + 0.5* cross3d((Rp - RprevVert),S/modS); % 0.5/N_v  * cross(S/modS, (RcurrVert - RprevVert)) +

% Loop over all vetices and add contribution due to presence of r_p in the area calculation of each triangular area
for i = 1:N_v
	v = vertices(i);
	if i < N_v
		vnext = vertices(i+1);
	else
		vnext = vertices(1);	
	end

	Rv = [vertexcoords(v,:),0];
	Rvnext = [vertexcoords(vnext,:),0];

	S_v = 0.5 * cross3d((Rp - Rvnext), (Rv - Rvnext));
	modS_v = norm(S_v);

	dA = dA + 0.5/N_v * cross3d((Rv - Rvnext),S_v/modS_v);
end

dA2d = dA(1:2);