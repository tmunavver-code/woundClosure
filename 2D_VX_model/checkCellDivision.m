function [celldata, param, division_happened] = checkCellDivision(celldata, param)
% CHECKCELLDIVISION Force-based cell division triggered when cell area exceeds
% Area_div_thresh, gated by Bell's law on area-elastic pressure and Metropolis energy.

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

    % Per-cell ka including wound margin stiffening
    ka_cell = param.ka;
    isWoundCell = isfield(param, 'cellIDtoContract') && ismember(cellID, param.cellIDtoContract);
    if isWoundCell
        if isfield(celldata, 'ka_wound_factor_current')
            ka_cell = ka_cell * celldata.ka_wound_factor_current;
        elseif isfield(param, 'ka_wound_factor')
            ka_cell = ka_cell * param.ka_wound_factor;
        end
    end

    % Area-elastic pressure driving force
    F_drive = 2.0 * ka_cell * (Acell - 1.0);

    % Effective force scale modulation
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
                % Metropolis uphill acceptance
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
                break; % One division per call
            else
                fprintf(1, 'DIV_REJECTED_DE: %f\n', dE);
            end
        end
    end
end

end
