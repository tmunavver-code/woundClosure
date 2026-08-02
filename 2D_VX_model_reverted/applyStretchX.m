
%% Applying stretching
deltaStretch = floor(0.5*param.Nsteps);

if param.ApplyStretchX == 1  && tstep - param.StretchAtStep <= deltaStretch
    param.Lx = param.Lx0;
    param.Ly = param.Ly0;
    if tstep >= param.StretchAtStep
        SR = 1 + (param.StretchRatio - 1)*(tstep - param.StretchAtStep)/deltaStretch; % stretching over 20 timesteps
        param.Lx = SR*param.Lx0;
        
        if param.BoxIncompressibility == 1
            param.Ly = param.Ly0/SR;
        end
        
    end
end