function f = getVertexForcesClassical(celldata,param,tstep)
%% getVertexForces: get vertex forces from gradient of an energy functional


f = zeros(size(celldata.f));

for cellID = 1:celldata.nCells
    
    Acell          = celldata.A(cellID);
    Pcell          = celldata.P(cellID);
    p0             = param.p0;
    rstiff         = param.rstiff;
    ka             = param.ka;
    vertexcoordsold   = celldata.r;
    Lx                = param.Lx;
    Ly                = param.Ly;
    vertices       = celldata.connec{cellID};
    fcell          = zeros(size(celldata.f));
    
    for i = 1:length(param.cellIDtoContract)
        if cellID == param.cellIDtoContract(i) && tstep > 55
            p0 = p0*param.multFactorForContraction;
        end
    end
    
    
    % Modify verecoords acocuntin for peridicity
    vertexcoords = modifyVerticesForPeriodicity(vertexcoordsold,vertices,Lx,Ly);
    
    
    for j = 1:length(vertices) % In anticlockwise order
        currVert = vertices(j);
        
        if j == length(vertices)
            nextVert = vertices(1);
        else
            nextVert = vertices(j+1);
        end
        
        if j == 1
            prevVert = vertices(end);
        else
            prevVert = vertices(j-1);
        end
        
        % Derivative of the cell perimeter with respect to the vertices (2 by 1)
        dPeri = getPerimeterDerivative(vertexcoords,currVert,nextVert,prevVert);
        
        % Derivative of the area with next vertex and the pericenter (2 by 1)
        dA    = getAreaDerivative(vertexcoords,vertices,currVert,nextVert,prevVert);
        
        if ismember(cellID,param.SelfPropellingCellIDs)
            Theta =  param.meanPropulsionAngle(1) + randn(1);
            polarityVector = [cos(Theta),sin(Theta)];
            % distribute the force among all vertices (the greater the
            % number of vertices, the lower the force per vertices)
            selfPropulsionForce = param.vel0 * polarityVector/length(celldata.connec{cellID});
        else
            selfPropulsionForce = [0, 0];
        end
        
        
        % Adding contributions from a cell to its vertices
        fcell(currVert,:) = fcell(currVert,:) - ...
            ka *2.0 * (Acell - 1) * dA - 2.0/rstiff * (Pcell - p0) * dPeri + selfPropulsionForce;
        
        
    end
    
    % Implementing fixed boundary
    if param.isBoundaryFixed == 1
       fcell(celldata.boundaryNodes,:) = 0; 
    end
    
    % Summing over all cells
    f = f + fcell;
end
