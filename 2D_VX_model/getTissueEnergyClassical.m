function totalenergy = getTissueEnergyClassical(celldata,param)
%% getTissueEnergyClassical: conservative energy used by topology gates

% get the total energy
totalenergy = 0;
for cellID = 1:celldata.nCells
    
    Acell  = celldata.A(cellID);
    Pcell  = celldata.P(cellID);
    p0     = param.p0;
    rstiff = param.rstiff;
    ka     = param.ka;
    
    isWoundCell = isfield(param, 'cellIDtoContract') && ...
        ismember(cellID, param.cellIDtoContract);
    if isWoundCell
        if isfield(celldata, 'ka_wound_factor_current')
            ka = ka * celldata.ka_wound_factor_current;
        elseif isfield(param, 'ka_wound_factor')
            ka = ka * param.ka_wound_factor;
        end
        if isfield(celldata, 'contractility_wound_current')
            rstiff = rstiff / celldata.contractility_wound_current;
        elseif isfield(param, 'contractility_wound')
            rstiff = rstiff / param.contractility_wound;
        end
    end
    
    % Cell energy
    energycell = ka*(Acell - 1.0)^2 + 1/rstiff * (Pcell - p0)^2;
    
    % Summing over all cells
    totalenergy = totalenergy + energycell;
end

% Wound-boundary cell energy (constant between topology changes, no vertex force).
if isfield(param, 'woundBoundaryCellEnergy') && ...
        param.woundBoundaryCellEnergy > 0 && ...
        isfield(param, 'cellIDtoContract')
    totalenergy = totalenergy + param.woundBoundaryCellEnergy * ...
        numel(unique(param.cellIDtoContract));
end

% Purse-string line energy (its gradient is the purse-string force).
if isfield(celldata, 'woundEdges') && ~isempty(celldata.woundEdges)
    if isfield(celldata, 'lambda_current')
        lambda = celldata.lambda_current;
    elseif isfield(param, 'lambda_purse_string')
        lambda = param.lambda_purse_string;
    else
        lambda = 0;
    end
    woundPerimeter = 0;
    for edgeID = 1:size(celldata.woundEdges, 1)
        d = celldata.r(celldata.woundEdges(edgeID,2),:) - ...
            celldata.r(celldata.woundEdges(edgeID,1),:);
        d(1) = d(1) - param.Lx * round(d(1) / param.Lx);
        d(2) = d(2) - param.Ly * round(d(2) / param.Ly);
        woundPerimeter = woundPerimeter + norm(d);
    end
    totalenergy = totalenergy + lambda * woundPerimeter;
end

% Optional wound-hole area elasticity (same polygon as the force term).
if isfield(param, 'enableWoundResistance') && param.enableWoundResistance && ...
        isfield(celldata, 'woundArea0') && celldata.woundArea0 > 0 && ...
        isfield(celldata, 'Kw_current') && celldata.Kw_current > 1e-6 && ...
        isfield(celldata, 'woundEdges') && ~isempty(celldata.woundEdges)
    [woundLoop, loopValid] = getOrderedWoundLoop(celldata.woundEdges);
    if loopValid
        woundData = celldata;
        woundData.connec{celldata.nCells + 1} = woundLoop;
        woundArea = getPolygonalCellArea(woundData, celldata.nCells + 1, param);
        totalenergy = totalenergy + celldata.Kw_current * ...
            (woundArea - celldata.woundArea0)^2;
    end
end

end
