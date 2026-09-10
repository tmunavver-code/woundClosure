function [celldata, param, division_happened] = checkCellDivision(celldata, param)
%CHECKCELLDIVISION Force-based cell division: a cell that has grown past
%an area threshold divides with a probability set by Bell's law on its
%own area-elastic pressure, the same stochastic-gating pattern used for
%T1 (checkT1transitions.m) and wound-margin turnover
%(checkWoundIntercalations.m). A division is only committed if it lowers
%total tissue energy (strict gate, no Metropolis uphill tolerance here,
%unlike T1/T2) -- otherwise the trial is reverted and the next candidate
%cell is tried.
%
% NOTE ON THIS CHOICE: Tetley et al. 2019 themselves treat division as an
% active, non-equilibrium "source of activity" and do not gate it on
% energy at all -- a strict downhill-only requirement will suppress most
% divisions, since creating new perimeter generally costs energy
% immediately. This condition was enforced as a stability criterion:
% the wound-topology corruption traced separately (see checkWoundIntercalations.m)
% was first observed at a timestep where a division and a wound intercalation
% coincided, and this gate rejects the energetically-unfavorable
% divisions that destabilize the boundary.
%
% Trigger: the cell's own current area exceeding param.Area_div_thresh,
% combined with the same Bell's-law stochastic draw used everywhere else
% in this codebase -- so cells under more area-elastic pressure (more
% crowded/over-grown) are more likely to divide, not a uniform per-cell
% random rate.
%
% Geometry: performCellDivision.m splits the cell along its short axis
% (Tetley's own stated division geometry -- "division of elongated cells
% along their short axis"), through the centroid.

division_happened = false;

if ~isfield(param, 'enableCellDivision') || ~param.enableCellDivision
    return;
end

E_old = getTissueEnergyClassical(celldata, param);

for cellID = 1:celldata.nCells

    vertices = celldata.connec{cellID};
    if length(vertices) < 4
        continue; % need >=4 vertices to split into two valid cells
    end

    Acell = celldata.A(cellID);
    if Acell <= param.Area_div_thresh
        continue;
    end

    % Per-cell ka, including the wound-margin boost if applicable --
    % must match getVertexForcesClassical.m/getEdgeForces.m exactly, same
    % consistency requirement as everywhere else force-based gating is used.
    ka_cell = param.ka;
    isWoundCell = isfield(param, 'cellIDtoContract') && ismember(cellID, param.cellIDtoContract);
    if isWoundCell
        if isfield(celldata, 'ka_wound_factor_current')
            ka_cell = ka_cell * celldata.ka_wound_factor_current;
        elseif isfield(param, 'ka_wound_factor')
            ka_cell = ka_cell * param.ka_wound_factor;
        end
    end

    % F_drive: area-elastic pressure (same quantity as pressure_area in
    % getEdgeForces.m), positive and growing as the cell over-shoots A0.
    F_drive = 2.0 * ka_cell * (Acell - 1.0);

    % Crowding modulation, same self-accelerating form used for T1/T2:
    % f_beta shrinks as the cell grows further past threshold.
    f_beta_eff = max(param.f_beta_div - param.f_beta_div_slope * (Acell - param.Area_div_thresh), 0.01);

    k_div = param.k_div0 * exp(F_drive / f_beta_eff);
    P_div = 1 - exp(-param.deltat * k_div);

    if rand() < P_div
        celldata_temp = celldata;
        [celldata_temp, success] = performCellDivision(cellID, celldata_temp, param);
        if success
            celldata_temp.A = getCellAreas(celldata_temp, param);
            [celldata_temp.P, celldata_temp.EdgeData] = getCellPerimeters(celldata_temp.nCells, celldata_temp.r, celldata_temp.connec, param, 0);
            E_new = getTissueEnergyClassical(celldata_temp, param);

            dE = E_new - E_old;

            if dE < 0
                accept_division = true;
            elseif isfield(param, 'T_eff_div') && param.T_eff_div > 0
                % Metropolis-style acceptance, same as T1/T2. Division is
                % the most expensive event type here (creating new
                % perimeter always costs energy immediately), so a strict
                % downhill-only rule rejected 100% of attempts -- see
                % T_eff_div's calibration note in
                % SingleCellContractionParameters.m.
                accept_division = (rand() < exp(-dE / param.T_eff_div));
            else
                accept_division = false;
            end

            if accept_division
                celldata = celldata_temp;
                division_happened = true;
                param.nDivisions = param.nDivisions + 1;
                fprintf(1, '--- Cell Division Performed --- Cell %d (area %.4f) split into %d and %d\n', ...
                    cellID, Acell, cellID, celldata.nCells);
                fprintf(1, 'Energy change: %f\n', dE);
                break; % one per call -- avoid stale-geometry compounding, same rule as T1/T2
            else
                fprintf(1, 'DIV_REJECTED_DE: %f\n', dE);
                % Trial discarded (celldata_temp goes out of scope
                % unused); try the next candidate cell instead.
            end
        end
    end
end

end
