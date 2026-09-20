function f = getVertexForcesClassical(celldata,param,tstep)
%% getVertexForcesClassical: Computes net forces on vertices from energy gradients
% Includes stability guards for degenerate cells (<3 vertices).

f = zeros(size(celldata.f));
Lx = param.Lx;
Ly = param.Ly;

% --- Active lamellipodial wound crawling force ---
haveWoundCrawl = isfield(param, 'enableWoundCrawling') && param.enableWoundCrawling && ...
    isfield(param, 'crawl_force0') && param.crawl_force0 > 0 && ...
    isfield(celldata, 'woundEdges') && ~isempty(celldata.woundEdges);
if haveWoundCrawl
    woundVerts = unique(celldata.woundEdges(:));
    refPt = celldata.r(woundVerts(1), :);
    dispFromRef = celldata.r(woundVerts, :) - refPt;
    dispFromRef(:,1) = dispFromRef(:,1) - Lx * round(dispFromRef(:,1) / Lx);
    dispFromRef(:,2) = dispFromRef(:,2) - Ly * round(dispFromRef(:,2) / Ly);
    woundCenter = refPt + mean(dispFromRef, 1);
    % Crawling force decays as purse-string recruits
    if isfield(celldata, 'crawl_current')
        current_crawl = celldata.crawl_current;
    else
        current_crawl = param.crawl_force0 * exp(-tstep / param.tau_decay_crawl);
    end
else
    current_crawl = 0;
end

for cellID = 1:celldata.nCells
    
    vertices = celldata.connec{cellID};
    
    % Guard clause: skip invalid cells (<3 vertices)
    if length(vertices) < 3
        continue;
    end
    
    Acell          = celldata.A(cellID);
    Pcell          = celldata.P(cellID);
    p0             = param.p0;
    rstiff         = param.rstiff;
    ka             = param.ka;
    vertexcoordsold   = celldata.r;
    fcell          = zeros(size(celldata.f));
    
    % Wound-cell specific contractility ramp
    isWoundCell = isfield(param, 'cellIDtoContract') && ismember(cellID, param.cellIDtoContract);
    if isWoundCell
        if isfield(celldata, 'ka_wound_factor_current')
            ka = ka * celldata.ka_wound_factor_current;
        elseif isfield(param, 'ka_wound_factor')
            ka = ka * param.ka_wound_factor;
        end
        if isfield(celldata, 'contractility_wound_current')
            rstiff = rstiff / celldata.contractility_wound_current;
        elseif isfield(param, 'contractility_wound')
            rstiff = rstiff / param.contractility_wound;  % lower rstiff = higher tension
        end
    end
    
    % Modify vertex coordinates accounting for periodicity
    vertexcoords = modifyVerticesForPeriodicity(vertexcoordsold,vertices,Lx,Ly);

    % Crawling force direction for margin cells
    if isWoundCell && current_crawl > 0
        cellBary = mean(vertexcoords, 1);
        dToWound = woundCenter - cellBary;
        dToWound(1) = dToWound(1) - Lx * round(dToWound(1) / Lx);
        dToWound(2) = dToWound(2) - Ly * round(dToWound(2) / Ly);
        distToWound = norm(dToWound);
        if distToWound > 1e-10
            crawlForce = current_crawl * (dToWound / distToWound) / length(vertices);
        else
            crawlForce = [0, 0];
        end
    else
        crawlForce = [0, 0];
    end

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
            ka *2.0 * (Acell - 1) * dA - 2.0/rstiff * (Pcell - p0) * dPeri + selfPropulsionForce + crawlForce;
        
    end
    
    if param.isBoundaryFixed == 1
       fcell(celldata.boundaryNodes,:) = 0; 
    end
    
    f = f + fcell;
end

% --- Edge-Specific Purse-String force ---
if isfield(celldata, 'woundEdges') && ~isempty(celldata.woundEdges)

    F_purse_string = zeros(size(f));
    % Recruit purse-string tension via feedback or decay ramp
    if isfield(celldata, 'lambda_current')
        current_lambda = celldata.lambda_current;
    else
        current_lambda = param.lambda_purse_string * min(1.0, tstep / param.t_ramp_purse_string);
    end

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

% --- Wound-hole elastic resistance (Tetley Kw) ---
if isfield(celldata, 'woundArea0') && celldata.woundArea0 > 0 && ...
        isfield(celldata, 'Kw_current') && celldata.Kw_current > 1e-6 && ...
        isfield(celldata, 'woundEdges') && ~isempty(celldata.woundEdges)

    [woundLoopVerts, woundLoopValid] = getOrderedWoundLoop(celldata.woundEdges);
    if woundLoopValid
        celldata_wound = celldata;
        pseudoID = celldata.nCells + 1;
        celldata_wound.connec{pseudoID} = woundLoopVerts;
        woundArea_current = getPolygonalCellArea(celldata_wound, pseudoID, param);

        woundVertCoords = modifyVerticesForPeriodicity(celldata.r, woundLoopVerts, Lx, Ly);
        nWV = length(woundLoopVerts);
        for j = 1:nWV
            currV = woundLoopVerts(j);
            if j == nWV
                nextV = woundLoopVerts(1);
            else
                nextV = woundLoopVerts(j+1);
            end
            if j == 1
                prevV = woundLoopVerts(end);
            else
                prevV = woundLoopVerts(j-1);
            end
            dA_wound = getAreaDerivative(woundVertCoords, woundLoopVerts, currV, nextV, prevV);
            f(currV,:) = f(currV,:) - celldata.Kw_current * 2.0 * (woundArea_current - celldata.woundArea0) * dA_wound;
        end
    end
end
% --- End wound-hole elastic resistance ---

end