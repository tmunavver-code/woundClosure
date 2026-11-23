function [celldata, stepScale] = updateVertexPositionsMonotonic(celldata, param)
% UPDATEVERTEXPOSITIONSMONOTONIC Backtracked overdamped vertex update.
% The explicit Euler step is halved until the mechanical energy does not increase.

energyBefore = getTissueEnergyClassical(celldata, param);
coordinatesBefore = celldata.r;
displacement = param.deltat * celldata.f / param.eta;
energyTolerance = 1e-12 * max(1, abs(energyBefore));

stepScale = 1;
for backtrack = 0:20
    trial = celldata;
    trial.r = coordinatesBefore + stepScale * displacement;
    trial.r(:,1) = mod(trial.r(:,1), param.Lx);
    trial.r(:,2) = mod(trial.r(:,2), param.Ly);
    trial.A = getCellAreas(trial, param);
    [trial.P, trial.EdgeData] = getCellPerimeters(trial.nCells, ...
        trial.r, trial.connec, param, 0);
    energyAfter = getTissueEnergyClassical(trial, param);

    if isfinite(energyAfter) && energyAfter <= energyBefore + energyTolerance
        celldata = trial;
        return;
    end
    stepScale = 0.5 * stepScale;
end

% No descent step found: keep the current configuration.
stepScale = 0;
celldata.r = coordinatesBefore;

end
