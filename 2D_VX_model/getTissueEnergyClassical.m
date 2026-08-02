function totalenergy = getTissueEnergyClassical(celldata,param)
%% getCurrentEnergy: returns total energyof the collection

% get the total energy
totalenergy = 0;
for cellID = 1:celldata.nCells
    
    Acell  = celldata.A(cellID);
    Pcell  = celldata.P(cellID);
    p0     = param.p0;
    rstiff = param.rstiff;
    ka     = param.ka;
    
    if cellID == param.cellIDtoContract
        rstiff = rstiff/param.multFactorForContraction;
    end
    
    % Energy pseudo-functionl
    energycell = ka*(Acell - 1.0)^2 + 1/rstiff * (Pcell - p0)^2;
    
    % Summing over all cells
    totalenergy = totalenergy + energycell;
end
