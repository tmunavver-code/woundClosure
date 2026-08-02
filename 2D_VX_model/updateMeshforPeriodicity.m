function [Vnew, Cnew] = updateMeshforPeriodicity(V, C, param)
% removes periodic images of vertices and updates the Cell connectivity
% information accordingly

Nvertices = length(V);
Ncells    = length(C);
Lx        = param.Lx;
Ly        = param.Ly;
TOL       = 1e-8;

masterlist = -ones(Nvertices,1);
Vperi      = V;
for vertID = 1:Nvertices
	Xpos = V(vertID,1);
	Ypos = V(vertID,2);
	perxflag = 0;
	peryflag = 0;
	
	% Recolate vertices outside of the box to their images inside the box
	if Xpos < 0
		Vperi(vertID,1) = Vperi(vertID,1) + Lx;
		perxflag = -1;
	elseif Xpos > Lx
		Vperi(vertID,1) = Vperi(vertID,1) - Lx;
		perxflag = 1;
	end
	
	if Ypos < 0
		Vperi(vertID,2) = Vperi(vertID,2) + Ly;
		peryflag = -1;
	elseif Ypos > Ly
		Vperi(vertID,2) = Vperi(vertID,2) - Ly;
		peryflag = 1;
	end
end

% Identify the vertices that are periodic images and assign one of them to
% be the master vertex that would be used as a degree of freedom
for vertID = 1:Nvertices
	Xpos = Vperi(vertID,1);
	Ypos = Vperi(vertID,2);
	peridiff = sqrt((Vperi(:,1) - Xpos).^2 + (Vperi(:,2) - Ypos).^2) < TOL;
	mastervertID = find(peridiff);
	if length(mastervertID) <= 1
		% Do nothing, self image identified
	else
		for j = 2:length(mastervertID)
			masterlist(mastervertID(j)) = mastervertID(1);
		end
	end
end

% Generate the new list of cells and vertices
Nverticesperi = length(find(masterlist==-1));
Vnew = zeros(Nverticesperi,2);
newvertID = 1;
oldnewmap = zeros(Nvertices,1);

for vertID = 1:Nvertices
	if masterlist(vertID) == -1
		Vnew(newvertID,:) = Vperi(vertID,:);
		oldnewmap(vertID) = newvertID;
		newvertID = newvertID + 1;
	end
end

% old new map of slave nodes based on masters
for vertID = 1:Nvertices
	if masterlist(vertID) ~= -1
		oldnewmap(vertID) = oldnewmap(masterlist(vertID));
	end
end

Cnew = C;
for cellID = 1:Ncells
	cellvertices = C{cellID};
	for j = 1:length(cellvertices)
		Cnew{cellID}(j) = oldnewmap(C{cellID}(j));
	end
end