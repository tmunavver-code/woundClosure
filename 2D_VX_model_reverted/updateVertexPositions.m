function celldata   = updateVertexPositions(celldata,param)
%% updateVertexPositions: update current vertex positions
% using the nodal force and viscosity coefficient, use backard Euler to
% update the nodal positions

f  = celldata.f;
r0 = celldata.r0;

deltat = param.deltat;
eta    = param.eta;

r = r0 + deltat*f/eta;

    
Lx = param.Lx;
Ly = param.Ly;
for vertID = 1:length(r)
    if r(vertID,1) > Lx
        r(vertID,1) = r(vertID,1) - Lx;
    end
    
    if r(vertID,1) < 0
        r(vertID,1) = r(vertID,1) + Lx;
    end
    
    if r(vertID,2) > Ly
        r(vertID,2) = r(vertID,2) - Ly;
    end
    
    if r(vertID,2) < 0
        r(vertID,2) = r(vertID,2) + Ly;
    end
end

celldata.r = r;