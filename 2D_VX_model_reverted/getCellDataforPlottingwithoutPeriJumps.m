function [Cnew, Vnew] = getCellDataforPlottingwithoutPeriJumps(nCells, Coords, Connec, param)

% Account for periodicity of cells when plotting, need to add vertices
% based on the periodicity information
V  = Coords;
C  = Connec;
Lx = param.Lx;
Ly = param.Ly;

Nverticesrep = 0;
for cellID = 1:nCells
    Nverticesrep = Nverticesrep + length(Connec{cellID});
end
Vnew = V;%zeros(Nverticesrep,2);

Cnew = C;
Nverticesold = length(V);
newvertexID  = Nverticesold;

n = 1;

for cellID = 1:nCells
    verticesmat = Connec{cellID};
    
    % As cells have vertices with a periodic jump, we take the first vertex
    % position and add new vertices on the same side as one of the vertices
    % (duplicating the vertices by adding their images)
    
    % taking the node closed to the center of mesh
    Vtemp = Vnew(verticesmat,:);
    dx = Vtemp(:,1) - param.Lx/2;
    dy = Vtemp(:,2) - param.Ly/2;
    disVec = sqrt(dx.^2 + dy.^2);
    [~,idx] = min(disVec);
    
    referenceVertex = verticesmat(idx);


    for j = 1:length(verticesmat)
        addnewVertflag = 0;
        Xnew = V(verticesmat(j),1);
        Ynew = V(verticesmat(j),2);
        
        if Xnew > V(referenceVertex,1) + Lx/2
            Xnew = Xnew - Lx;
            addnewVertflag = 1;
        end
        
        if Xnew < V(referenceVertex,1) - Lx/2
            Xnew = Xnew + Lx;
            addnewVertflag = 1;
        end
        
        if Ynew > V(referenceVertex,2) + Ly/2
            Ynew = Ynew - Ly;
            addnewVertflag = 1;
        end
        
        if Ynew < V(referenceVertex,2) - Ly/2
            Ynew = Ynew + Ly;
            addnewVertflag = 1;
        end
        
        
        % Add new vertex to the list
        if addnewVertflag == 1
             
            Vnew = [Vnew; Xnew, Ynew];
            newvertexID = newvertexID + 1;
            Cnew{cellID}(j) = newvertexID; % replaces the old vertices by the new ones
        end
        
        n = 1;
    end
    
    
end
