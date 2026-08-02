function Vnew = modifyVerticesForPeriodicity(V,verticesmat,Lx,Ly)
% update vertices coordinates based on peridicity for the currnt cell
% vertices

Vnew = V;
for j = 2:length(verticesmat)
	vertID = verticesmat(j);
	Xnew = V(vertID,1);
	Ynew = V(vertID,2);
	
	if Xnew > V(verticesmat(1),1) + Lx/2
		Xnew = Xnew - Lx;
	end
	
	if Xnew < V(verticesmat(1),1) - Lx/2
		Xnew = Xnew + Lx;
	end
	
	if Ynew > V(verticesmat(1),2) + Ly/2
		Ynew = Ynew - Ly;
	end
	
	if Ynew < V(verticesmat(1),2) - Ly/2
		Ynew = Ynew + Ly;
	end
	
	Vnew(vertID,1) = Xnew;
	Vnew(vertID,2) = Ynew;
end