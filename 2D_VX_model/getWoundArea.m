function woundArea = getWoundArea(celldata, param)
% GETWOUNDAREA Computes total wound area via periodic box area conservation:
% woundArea = Lx * Ly - sum(cell areas).

woundArea = param.Lx * param.Ly - sum(celldata.A);

if woundArea < 0
    woundArea = 0; % numerically closed
end

end
