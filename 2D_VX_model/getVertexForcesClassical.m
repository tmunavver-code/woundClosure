function f = getVertexForcesClassical(celldata,param,tstep)
%% getVertexForces: get vertex forces from gradient of an energy functional
%
% MODIFIED (Stable V2):
% - Added a guard clause to skip force calculations for invalid cells
%   (with < 3 vertices). This prevents the simulation from "blowing up".
%

f = zeros(size(celldata.f));
Lx = param.Lx;
Ly = param.Ly;

for cellID = 1:celldata.nCells
    
    vertices = celldata.connec{cellID};
    
    % --- NEW GUARD CLAUSE ---
    % If a 3-cell intercalation has reduced this cell to a 2-vertex
    % "line", it is no longer a valid cell. Skip all physics for it.
    if length(vertices) < 3
        continue; % Skip to the next cellID
    end
    % --- END NEW GUARD CLAUSE ---
    
    Acell          = celldata.A(cellID);
    Pcell          = celldata.P(cellID);
    p0             = param.p0;
    rstiff         = param.rstiff;
    ka             = param.ka;
    vertexcoordsold   = celldata.r;
    fcell          = zeros(size(celldata.f));
    
    % (No p0 modification loop, this is correct)
    
    % Modify verecoords acocuntin for peridicity
    vertexcoords = modifyVerticesForPeriodicity(vertexcoordsold,vertices,Lx,Ly);
    
    
    for j = 1:length(vertices) % In anticlockwise order
        % ... (rest of the function is identical) ...
        
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
        
        dPeri = getPerimeterDerivative(vertexcoords,currVert,nextVert,prevVert);
        dA    = getAreaDerivative(vertexcoords,vertices,currVert,nextVert,prevVert);
        
        if ismember(cellID,param.SelfPropellingCellIDs)
            Theta =  param.meanPropulsionAngle(1) + randn(1);
            polarityVector = [cos(Theta),sin(Theta)];
            selfPropulsionForce = param.vel0 * polarityVector/length(celldata.connec{cellID});
        else
            selfPropulsionForce = [0, 0];
        end
        
        fcell(currVert,:) = fcell(currVert,:) - ...
            ka *2.0 * (Acell - 1) * dA - 2.0/rstiff * (Pcell - p0) * dPeri + selfPropulsionForce;
        
    end
    
    if param.isBoundaryFixed == 1
       fcell(celldata.boundaryNodes,:) = 0; 
    end
    
    f = f + fcell;
end

% ... (The NEW Edge-Specific Purse-String force section remains unchanged) ...
if isfield(celldata, 'woundEdges') && ~isempty(celldata.woundEdges)
    
    F_purse_string = zeros(size(f));
    current_lambda = param.lambda_purse_string * min(1.0, tstep / param.t_ramp_purse_string);

    for i = 1:size(celldata.woundEdges, 1)
        v1_idx = celldata.woundEdges(i, 1);
        v2_idx = celldata.woundEdges(i, 2);
        r1 = celldata.r(v1_idx, :);
        r2 = celldata.r(v2_idx, :);
        
        vec_1_to_2 = r2 - r1;
        vec_1_to_2(1) = vec_1_to_2(1) - Lx * round(vec_1_to_2(1) / Lx);
        vec_1_to_2(2) = vec_1_to_2(2) - Ly * round(vec_1_to_2(2) / Ly);
        
        dist = norm(vec_1_to_2);
        
        if dist > 1e-6 
            unit_vec = vec_1_to_2 / dist;
            force_on_v1 = current_lambda * unit_vec;
            force_on_v2 = -current_lambda * unit_vec;
            F_purse_string(v1_idx, :) = F_purse_string(v1_idx, :) + force_on_v1;
            F_purse_string(v2_idx, :) = F_purse_string(v2_idx, :) + force_on_v2;
        end
    end
    
    f = f + F_purse_string;
    
end

end