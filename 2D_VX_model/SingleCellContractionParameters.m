% T1 transition parameters
param.maxT1relaxsteps = 30;
param.nT1      = 0;         % to store the number of T1 transitions occured

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
% A wound contact is also eligible when it is short relative to the
% typical edge length of its owning cell. Handoffs are tried before trims.
param.L_intercalation_ratio_thresh = 0.75;
% Start a wound-margin T1 from an almost collapsed edge.  The replacement
% edge then grows under the forces instead of being inserted at a finite
% length, avoiding an artificial energy jump at the topology change.
param.woundT1NewEdgeLength = 0.01;
% The T2 Bell gate relaxes as the wound perimeter contracts.  The upper
% value suppresses weak early events; the lower value is just below the
% zero-force probability 1-exp(-k_off0*deltat), so geometry can trigger
% late handoffs after the purse string has shortened the margin.
param.bellProbabilityThresholdT2 = 0.0125;
param.bellProbabilityThresholdT2Min = 0.0090;
param.minWoundBoundaryCells = 5;
% Discrete free-energy cost for maintaining one cell on the wound margin.
% This term has no positional derivative and therefore adds no resistance
% force; it supplies an energetic drive for gradual boundary-cell exit.
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
% Let boundary handoffs reduce the margin before opposing vertices merge.
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
% Deterministic Bell gate. At the default dt and base rate, the zero-force
% probability is 1-exp(-0.01)=0.00995. A 0.02 cutoff suppresses baseline
% events while retaining transitions whose probability is raised by force.
param.bellProbabilityThreshold = 0.02;
% Bell's law characteristic force scales
param.f_beta_T1 = 22.271433;
param.f_beta_T2 = 0.494935;

% Metropolis uphill acceptance tolerances
param.T_eff_T1 = 0;
param.T_eff_T2 = 0;
param.f_beta_length_slope_T1 = 1.0;
param.f_beta_length_slope_T2 = 1.0;

% plotting parameters
param.nVisualisation = 50;     % frequency of current state visualisation

% simulation time parameters
param.deltat   = 0.01;      % initial timestep
param.Tsim     = 0;         % simulation time
param.Nsteps   = 3000;       % Increased for wound closure (was 1000)
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


% --- Edge-Specific Purse-String Parameters ---
param.lambda_purse_string = 4.0;  % Increased for wound closure (was 2.0)
param.t_ramp_purse_string = 50;   % Faster ramp-up (was 100)

% --- Wound-Cell Contractility Parameters (margin cells only) ---
param.ka_wound_factor     = 2.0;  % wound cells have 2x area stiffness
param.contractility_wound = 3.0;  % wound cells have 3x perimeter contractility

% --- Wound-Directed Crawling Parameters ---
param.enableWoundCrawling = true;
param.crawl_force0    = 3.0;
param.tau_decay_crawl = 100;

% --- Force-Based Closure Strategy Parameters ---
param.Kw0    = param.ka * 2;   % initial wound-hole area-elastic modulus
param.tau_Kw = 2;              % decay time for Kw
param.k_recruit = 1.0;         % rate constant for recruitment ODEs
% Kw remains available as a recruitment clock, but it does not exert a
% wound-area force or contribute a wound-area elastic-energy term.
param.enableWoundResistance = false;

% With fixed coefficients and conservative forces, the overdamped vertex
% dynamics should relax down one mechanical-energy landscape. Topology
% changes are therefore restricted to energy-lowering events above.
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
