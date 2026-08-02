% T1 transition parameters
param.maxT1relaxsteps = 30;
param.nT1      = 0;         % to store the number of T1 transitions occured

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%% NEW: 3-Cell Intercalation Parameters %%%%%%%%%%%%%%%%%%%%
% Implements the T2-like extrusion/intercalation at the wound margin.
% This is a separate event from T1 (4-cell) transitions.
param.enableWoundIntercalations = true;
param.L_intercalation_thresh = 0.5; % Threshold for wound-edge shrink
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%%%%%%%%%%%%%%% T1 Intercalation (4-Cell) Parameters %%%%%%%%%%%%%%%%%%%%
param.enableT1transitions = false;   % MASTER TOGGLE to turn all T1 swaps on/off
param.enhanceT1atWound = false;    % Toggle for higher intercalation at wound edge
param.T1_TOL_bulk = 0.1;      % tolerance for bulk tissue
param.T1_TOL_wound = 0.2;     % permissive tolerance for wound-edge cells
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% plotting parameters
param.nVisualisation = 50;     % frequency of current state visualisation

% simulation time parameters
param.deltat   = 0.01;      % initial timestep
param.Tsim     = 0;         % simulation time
param.Nsteps   = 1000;       % simulation time (You may want to increase this)
param.isBoundaryFixed = 0;  % fixed boundary doesn't allow for periodic jumps

% cells parameters
param.rstiff   = .5;               % stiffness factor
param.ka       = 1.0;               % area term premultiplier
param.eta      = 1;               % vertex viscosity

% ... (other parameters like self-propulsion, stretch, etc.) ...
param.vel0     = 0.0;
param.SelfPropellingCellIDs = 0;
param.StretchAtStep = 100000;
param.StretchRatio  = 1;
param.ApplyStretchX = 0;
param.BoxIncompressibility = 1;
param.Lx0 = param.Lx;
param.Ly0 = param.Ly;
param.cellIDstoTrack = 0;

%%%%%%%%%%%%%%%%%%%%%% Cell contraction parameter (OLD - DISABLED) %%%%%%%%%%%%%%%%%%
%     param.cellIDtoContract = [67 66 50 122 79 17 35 20 102 33 31 144 ...
%         57 59 63 26 104 123 94 64 24 77 32 59 ...
%         57 63 24 81 130 99 37 90 128];
%     param.multFactorForContraction = .25;
%%%%%%%%%%%%%%%%%%%%%% Cell contraction parameter (OLD - DISABLED) %%%%%%%%%%%%%%%%%%


% %%%%%%%%%%%%%%%%%%%%% NEW: Edge-Specific Purse-String Parameters %%%%%%%%%%%%%%%%%%
param.lambda_purse_string = 2.0;  % Max tension of the purse-string cable
param.t_ramp_purse_string = 100;  % Number of timesteps to ramp up to max tension
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%%%%%%%%%%%%%%%%%%% NEW: Force Visualization Parameters %%%%%%%%%%%%%%%%%%%%%%%%%
param.plotForces = true;      % Set to true to enable force plotting
param.plotForces_Interval = 20; % Plot forces every 20 timesteps
param.plotForces_FigHandle = 10; % Figure window to use for plotting
param.plotForces_Scale = 2.0;  % Tune this to make force arrows larger/smaller
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Target shape index
% ... (rest of the file) ...
if strcmp(Case,'fluid')
    fprintf(1,'Tissue is now %s\n',Case);
    param.p0 = 3.9;
    % parm.multFactorForContraction = min(param.multFactorForContraction,1.5); % No longer needed
elseif strcmp(Case,'solid')
    param.p0 = 2;
    % parm.multFactorForContraction = min(param.multFactorForContraction,2.5); % No longer needed
    fprintf(1,'Tissue is now %s\n',Case);
else
    error('Please choose a given tissue fluidity case : "fluid" or "solid".')
end