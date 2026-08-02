% T1 transition parameters
param.maxT1relaxsteps = 20;
param.T1_TOL   = 0.05;      % tolerance (edge length) to perform a T1 swap
param.nT1      = 0;         % to store the number of T1 transitions occured

% plotting parameters
param.nVisualisation = 50;     % frequency of current state visualisation

% simulation parameters
param.deltat   = 0.01;      % initial timestep
param.Tsim     = 0;         % simulation time
param.Nsteps   = 5000;      % simulation time
param.isBoundaryFixed = 0;  % fixed boundary doesn't allow for periodic jumps


% cell parameters
param.rstiff   = 0.5;               % stiffness factor
param.ka       = 1;               % area term premultiplier
param.eta      = 1;               % vertex viscosity

%%%%%%%%%%%%%%%%%%%%%% self-propulsion parameter %%%%%%%%%%%%%%%%%%%%%%

    param.SelfPropellingCellIDs = [12];
    param.vel0     = 8.0;               % velocity of self-propulsion
    param.meanPropulsionAngle = [pi/2,3*pi/2];

%%%%%%%%%%%%%%%%%%%%%% self-propulsion parameter %%%%%%%%%%%%%%%%%%%%%%

% stretchting parameters
param.StretchAtStep = 100000;
param.StretchRatio  = 1;
param.ApplyStretchX = 0;       % apply stretch in the x direction
param.BoxIncompressibility = 1;      % conserve area when applying stretch
param.Lx0 = param.Lx;
param.Ly0 = param.Ly;
param.cellIDstoTrack = 0;   % tracking neighbouring cells

% Cell contraction parameter
param.cellIDtoContract = 0;

% Target shape index
% threhsold between solid and fluid phase is at ~ 3.8
if strcmp(Case,'fluid')
    fprintf(1,'Tissue is now %s\n',Case);
    param.p0 = 3.9; % below 4 for better simulation results
elseif strcmp(Case,'solid')
    param.p0 = 3;
    fprintf(1,'Tissue is now %s\n',Case);
else
    error('Please choose a given tissue fluidity case : "fluid" or "solid".')
end