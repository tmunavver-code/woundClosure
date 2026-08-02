% T1 transition parameters
param.maxT1relaxsteps = 30;
param.T1_TOL   = 0.1;      % tolerance to perform a T1 swap
param.nT1      = 0;        % number of T1 transitions

% plotting parameters
param.nVisualisation = 20;  % frequency of current state visualisation

% simulation parameters
param.deltat   = 0.01;      % initial timestep
param.Tsim     = 0;         % simulation time
param.Nsteps   = 1000;      % simulation time
param.isBoundaryFixed = 0;  % fixed boundary doesn't allow for periodic jumps

% cell parameters
param.rstiff   = 0.5;        % stiffness factor
param.ka       = 1.0;        % area term premultiplier
param.eta      = 1.0;        % vertex viscosity

% self-propulsion parameter
param.vel0     = 0.0;               % velocity of self-propulsion
param.SelfPropellingCellIDs = 0;

%%%%%%%%%%%%%%%%%%%%%% stretchting parameters %%%%%%%%%%%%%%%%%%%%%%

    param.StretchAtStep = floor(param.Nsteps/5);
    param.StretchRatio  = 2;
    param.ApplyStretchX = 1;             % apply stretch in the x direction
    param.BoxIncompressibility = 1;      % conserve area when applying stretch
    param.Lx0 = param.Lx;
    param.Ly0 = param.Ly;
    param.cellIDstoTrack = [13,33,12,28];   % tracking neighbouring cells

%%%%%%%%%%%%%%%%%%%%%% stretchting parameters %%%%%%%%%%%%%%%%%%%%%%

% Cell contraction parameter
param.cellIDtoContract = 0;

% Target shape index
% threhsold between solid and fluid phase is at ~ 3.8
if strcmp(Case,'fluid')
    fprintf(1,'Tissue is now %s\n',Case);
    param.p0 = 3.9; % below 4 for better simulation results
elseif strcmp(Case,'solid')
    param.p0 = 3.6;
    fprintf(1,'Tissue is now %s\n',Case);
else
    error('Please choose a given tissue fluidity case : "fluid" or "solid".')
end