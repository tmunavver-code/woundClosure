% Default model parameters for the wound-closure simulations.
% Requires param.Lx, param.Ly and Case ('fluid' or 'solid') to be set first.

% T1 transition parameters
param.maxT1relaxsteps = 30;
param.nT1      = 0;         % counter of accepted T1 transitions

% --- Force-Based Cell Division Parameters ---
param.enableCellDivision = false;
param.allowWoundMarginDivision = false;
param.nDivisions = 0;
param.Area_div_thresh = 1.102315;
param.k_div0 = 1.0;            % Base division rate
param.f_beta_div = 0.250021;
param.f_beta_div_slope = 1.0;
% Metropolis tolerance for division energy gate
param.T_eff_div = 3.512995;

% --- Wound Intercalation (3-Cell) Parameters ---
param.enableWoundIntercalations = true;
param.L_intercalation_thresh = 0.5; % Threshold for wound-edge shrink
% A wound contact is also eligible when short relative to its cell's mean edge length.
param.L_intercalation_ratio_thresh = 0.75;
% Length of the new edge created by a wound-margin T1 (grows under the forces).
param.woundT1NewEdgeLength = 0.01;
% Bell threshold for margin events; relaxes from the upper to the lower value
% as the wound perimeter contracts.
param.bellProbabilityThresholdT2 = 0.0125;
param.bellProbabilityThresholdT2Min = 0.0090;
param.minWoundBoundaryCells = 5;
% Energy per wound-boundary cell (no positional force; favours margin turnover).
param.woundBoundaryCellEnergy = 5.0;
param.nWoundIntercalations = 0;
param.nWoundCornerTrims = 0;

% --- Wound Sealing Parameters ---
param.enableWoundSealing = false;
param.nSeals = 0;
param.L_seal_thresh = 0.15;
param.k_seal0 = 1.0;
param.f_beta_seal = 2.255200;
param.f_beta_seal_slope = 1.0;
param.bellProbabilityThresholdSeal = 0.02;
% Sealing starts once the margin has at most this many cells.
param.maxBoundaryCellsForSealing = 5;
param.maxSealVertexDegree = 5;
% Metropolis tolerance for wound sealing energy gate
param.T_eff_seal = 0;

% --- T1 Intercalation (4-Cell) Parameters ---
param.enableT1transitions = true;    % Master toggle for T1 transitions
param.enhanceT1atWound = true;       % Toggle for higher intercalation at wound edge
param.T1_TOL_bulk = 0.05;            % Tolerance for bulk tissue
param.T1_TOL_wound = 0.15;           % Tolerance for wound-edge cells


% --- Force-Based Transition Parameters ---
param.useForceBasedTransitions = true;
param.k_off0 = 1.0;
% Deterministic Bell threshold (zero-force probability is 1-exp(-0.01) = 0.00995).
param.bellProbabilityThreshold = 0.02;
% Bell's law characteristic force scales
param.f_beta_T1 = 22.271433;
param.f_beta_T2 = 0.494935;

% Metropolis temperatures (0 = accept only energy-lowering events)
param.T_eff_T1 = 0;
param.T_eff_T2 = 0;
param.f_beta_length_slope_T1 = 1.0;
param.f_beta_length_slope_T2 = 1.0;

% plotting parameters
param.nVisualisation = 50;     % frequency of current state visualisation

% simulation time parameters
param.deltat   = 0.01;      % initial timestep
param.Tsim     = 0;         % simulation time
param.Nsteps   = 3000;       % default length; drivers set their own
param.isBoundaryFixed = 0;  % fixed boundary doesn't allow for periodic jumps

% cells parameters
param.rstiff   = .5;               % stiffness factor
param.ka       = 1.0;               % area term premultiplier
param.eta      = 1;               % vertex viscosity

% Self-propulsion, stretching and tracking (inactive in the wound runs)
param.vel0     = 0.0;
param.SelfPropellingCellIDs = 0;
param.StretchAtStep = 100000;
param.StretchRatio  = 1;
param.ApplyStretchX = 0;
param.BoxIncompressibility = 1;
param.Lx0 = param.Lx;
param.Ly0 = param.Ly;
param.cellIDstoTrack = 0;



% --- Edge-Specific Purse-String Parameters ---
param.lambda_purse_string = 4.0;  % target purse-string line tension
param.t_ramp_purse_string = 50;   % ramp length (steps) when no recruitment clock is used

% --- Wound-Cell Contractility Parameters (margin cells only) ---
param.ka_wound_factor     = 2.0;  % wound cells have 2x area stiffness
param.contractility_wound = 3.0;  % wound cells have 3x perimeter contractility

% --- Wound-Directed Crawling Parameters ---
param.enableWoundCrawling = true;
param.crawl_force0    = 3.0;
param.tau_decay_crawl = 100;

% --- Recruitment clock and optional wound-area resistance ---
param.Kw0    = param.ka * 2;   % initial wound-hole area-elastic modulus
param.tau_Kw = 2;              % decay time for Kw
param.k_recruit = 1.0;         % rate constant for recruitment ODEs
% false: Kw only paces recruitment; true: Kw also acts as wound-area elasticity.
param.enableWoundResistance = false;

% true: fixed coefficients, no crawling, backtracking integrator.
param.enforceMonotonicEnergy = false;

% --- Force Visualization Parameters ---
param.plotForces = true;      % Set to true to enable force plotting
param.plotForces_Interval = 20; % Plot forces every 20 timesteps
param.plotForces_FigHandle = 10; % Figure window to use for plotting
param.plotForces_Scale = 2.0;  % Tune this to make force arrows larger/smaller


% Target shape index
if strcmp(Case,'fluid')
    fprintf(1,'Tissue is now %s\n',Case);
    % Target shape index p0 for fluid regime (Bi et al. 2015)
    param.p0 = 4.0;
elseif strcmp(Case,'solid')
    param.p0 = 2;
    fprintf(1,'Tissue is now %s\n',Case);
else
    error('Please choose a given tissue fluidity case : "fluid" or "solid".')
end
