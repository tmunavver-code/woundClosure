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

% --- Wound-directed crawling force (Trepat et al. 2014, Nat. Phys.) ---
% Early wound closure is driven by outward-pointing lamellipodial
% protrusion of margin cells into the wound (the OPTL), which decays as
% the purse-string ring (IPTL) takes over. This is a separate, active
% force -- distinct from the passive purse-string line tension below --
% oriented from each margin cell's barycenter toward the wound centroid.
% See README_force_based_transitions.md for the rationale.
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
    % FORCE-BASED STRATEGY (Phase F1): crawl force now decays as purse-string
    % tension actually builds (celldata.crawl_current, integrated in
    % myMainClassical.m), not on an independent fixed clock. Falls back to
    % the old fixed-decay formula only if that state isn't present (e.g.
    % older test scripts), so nothing breaks that isn't updated.
    if isfield(celldata, 'crawl_current')
        current_crawl = celldata.crawl_current;
    else
        current_crawl = param.crawl_force0 * exp(-tstep / param.tau_decay_crawl);
    end
else
    current_crawl = 0;
end
% --- End wound-directed crawling force setup ---

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
    
    % --- Wound-cell specific contractility ---
    % FORCE-BASED STRATEGY (Phase F5): the boost now ramps in via
    % celldata.ka_wound_factor_current / celldata.contractility_wound_current
    % (integrated in myMainClassical.m, gated on the wound's Kw having
    % decayed -- see force_based_closure_strategy.md), instead of applying
    % the full multiplier instantly at t=0. Falls back to the old instant
    % multiplier if that state isn't present.
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
    % --- End wound contractility ---
    
    % (No p0 modification loop, this is correct)
    
    % Modify verecoords acocuntin for peridicity
    vertexcoords = modifyVerticesForPeriodicity(vertexcoordsold,vertices,Lx,Ly);

    % --- Crawling force direction for this cell (if it is a wound-margin cell) ---
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
    % FORCE-BASED STRATEGY (Phase F1): tension now recruited via feedback
    % (celldata.lambda_current, integrated in myMainClassical.m, gated on
    % Kw decay) instead of a fixed tstep ramp. Falls back to the old ramp
    % if that state isn't present.
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

% --- Wound-hole elastic resistance (Phase F2, Tetley et al. 2019 Kw) ---
% The ablated hole retains a residual area-elastic resistance pulling it
% back toward its area right after ablation (celldata.woundArea0), with
% modulus celldata.Kw_current decaying to 0 over tau_Kw (integrated in
% myMainClassical.m). While Kw is large this resists deformation in
% *either* direction; combined with released bulk pre-stress this is what
% produces the initial passive expansion phase, before it fades and lets
% the purse-string close the wound. See force_based_closure_strategy.md.
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