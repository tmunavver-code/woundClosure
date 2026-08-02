function f = getVertexForcesClassical(celldata,param,tstep)
%% getVertexForces: get vertex forces from gradient of an energy functional
%
% MODIFIED:
% 1. REMOVED: Old whole-cell contraction (p0 modification).
% 2. ADDED: New edge-specific purse-string force calculation.
%

f = zeros(size(celldata.f));
Lx = param.Lx;
Ly = param.Ly;

for cellID = 1:celldata.nCells
    
    Acell          = celldata.A(cellID);
    Pcell          = celldata.P(cellID);
    p0             = param.p0;
    rstiff         = param.rstiff;
    ka             = param.ka;
    vertexcoordsold   = celldata.r;
    vertices       = celldata.connec{cellID};
    fcell          = zeros(size(celldata.f));
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % --- OLD CONTRACTION MECHANISM (REMOVED) ---
    % for i = 1:length(param.cellIDtoContract)
    %     if cellID == param.cellIDtoContract(i) && tstep > 55
    %         p0 = p0*param.multFactorForContraction;
    %     end
    % end
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
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

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- NEW: PRIORITY 1 - EDGE-SPECIFIC PURSE-STRING FORCE ---
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This force is added *after* all the cell-based elastic forces are calculated.
% It acts *only* on the edges stored in celldata.woundEdges.

if isfield(celldata, 'woundEdges') && ~isempty(celldata.woundEdges)
    
    % Initialize a separate force matrix for the purse-string
    F_purse_string = zeros(size(f));
    
    % Get current tension, ramping up from t=0 to t=t_ramp
    current_lambda = param.lambda_purse_string * min(1.0, tstep / param.t_ramp_purse_string);

    % Loop over every edge segment in the wound margin
    for i = 1:size(celldata.woundEdges, 1)
        % Get vertex indices for this edge
        v1_idx = celldata.woundEdges(i, 1);
        v2_idx = celldata.woundEdges(i, 2);

        % Get vertex coordinates
        r1 = celldata.r(v1_idx, :);
        r2 = celldata.r(v2_idx, :);

        % Calculate the vector between them using Minimum Image Convention (MIC)
        % This correctly handles edges that cross the periodic boundary
        vec_1_to_2 = r2 - r1;
        vec_1_to_2(1) = vec_1_to_2(1) - Lx * round(vec_1_to_2(1) / Lx);
        vec_1_to_2(2) = vec_1_to_2(2) - Ly * round(vec_1_to_2(2) / Ly);
        
        dist = norm(vec_1_to_2);
        
        if dist > 1e-6 % Avoid division by zero
            unit_vec = vec_1_to_2 / dist;

            % Calculate the contractile force vector
            % Pulls v1 towards v2, and v2 towards v1
            force_on_v1 = current_lambda * unit_vec;
            force_on_v2 = -current_lambda * unit_vec;

            % Add this force to the purse-string force matrix
            F_purse_string(v1_idx, :) = F_purse_string(v1_idx, :) + force_on_v1;
            F_purse_string(v2_idx, :) = F_purse_string(v2_idx, :) + force_on_v2;
        end
    end
    
    % Add the purse-string forces to the total elastic forces
    f = f + F_purse_string;
    
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- END OF NEW SECTION ---
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end