% T1 transition parameters
param.maxT1relaxsteps = 30;
param.nT1      = 0;         % to store the number of T1 transitions occured

% %%%%%%%%%%%%%%%%%%%%% Force-Based Cell Division %%%%%%%%%%%%%%%%%%%%%%%%%%%
% Third topological event type, alongside T1 (checkT1transitions.m) and
% wound-margin turnover (checkWoundIntercalations.m). Tetley et al. 2019
% use cell division (age/fixed-rate) as one of two non-equilibrium
% "activity" sources; here the trigger is force-based instead: a cell
% whose area exceeds Area_div_thresh divides with a Bell's-law
% probability set by its own area-elastic pressure (checkCellDivision.m),
% no energy gate (division is meant to be an active, energy-costing
% process, not a passive relaxation -- see checkCellDivision.m's own
% comment). Geometry: split along the short axis through the centroid
% (performCellDivision.m), matching Tetley's stated division orientation.
% Calibration run: 1000 steps, fluid regime, division OFF, sampling every
% cell's area every 5 steps (n=42,600) to measure the tissue's own
% natural area distribution -- mean=1.031, median=1.040, p90=1.102,
% p95=1.125, p99=1.159, max=1.387. Area_div_thresh set at the 90th
% percentile (only the top ~10% largest cells at any moment are
% candidates); f_beta_div set at the median area-elastic pressure among
% cells already above that threshold (n=4260), same "characteristic
% force scale" logic as f_beta_T1/f_beta_T2.
param.enableCellDivision = true;
param.nDivisions = 0;
param.Area_div_thresh = 1.102315;
param.k_div0 = 1.0;            % no data-driven basis exists for this base rate, same limitation as k_off0
param.f_beta_div = 0.250021;
param.f_beta_div_slope = 1.0;
% Metropolis tolerance for division's energy gate, same pattern as
% T_eff_T1/T_eff_T2. Without it, a strict downhill-only rule rejected
% 807/807 attempts in a 1000-step run -- division existed but never
% fired. Calibrated from that run's own logged DIV_REJECTED_DE
% population: n=807, median=3.512995, mean=3.539929 (p10=2.67,
% p90=4.56). Note this is ~2x T_eff_T2 and ~90x T_eff_T1: creating new
% perimeter is genuinely the most energetically expensive topological
% event in this model, so it needs the largest tolerance to fire at all.
param.T_eff_div = 3.512995;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%%%%%%%%%%%%%%%% NEW: 3-Cell Intercalation Parameters %%%%%%%%%%%%%%%%%%%%
% Implements the T2-like extrusion/intercalation at the wound margin.
% This is a separate event from T1 (4-cell) transitions.
param.enableWoundIntercalations = true;
param.L_intercalation_thresh = 0.5; % Threshold for wound-edge shrink

% %%%%%%%%%%%%%%%%%%%%% Wound Sealing (checkWoundSealing.m) %%%%%%%%%%%%%%%%%
% Fourth topological event type. Merges margin vertices from OPPOSITE
% sides of the collapsing wound, converting free edges into ordinary
% shared edges -- this is what actually retires cells from the wound
% boundary at the end of closure. Without it the wound area reaches zero
% while the margin still reports its full cell count (measured: 0.00%
% area at step 700 with 15 margin cells / 22 free edges still present),
% which is the mismatch against the experimental videos.
% Gated by distance + Bell's law on the closing force + Metropolis energy
% acceptance, same as T1 / wound-margin turnover / division.
param.enableWoundSealing = true;
param.nSeals = 0;
% Calibrated by logging every cross-gap margin-vertex pair (no shared
% cell) over a 1000-step run, n=133,541. The closing force sharpens
% strongly as the gap narrows, which is what justifies gating on it:
%   gap<0.50: median F_drive=-0.05, 50% positive  (no net closing -- noise)
%   gap<0.20: median F_drive= 0.95, 54% positive
%   gap<0.15: median F_drive= 2.26, 66% positive
%   gap<0.10: median F_drive= 4.59, 85% positive
% L_seal_thresh=0.15 is where closing force clearly dominates while
% keeping a usable candidate pool; it also coincides with T1_TOL_wound,
% i.e. the same length scale at which junctions rearrange elsewhere.
% f_beta_seal = median closing force among those candidates (n=1876),
% same characteristic-force-scale rule as f_beta_T1/T2/div.
param.L_seal_thresh = 0.15;
param.k_seal0 = 1.0;
param.f_beta_seal = 2.255200;
param.f_beta_seal_slope = 1.0;
% T_eff_seal from the same run with the strict downhill-only gate, taking
% the median of the logged SEAL_REJECTED_DE population (n=218,
% median=0.090348, mean=0.137712). Sealing is by far the CHEAPEST
% topological event here -- ~0.09 versus 1.87 (T2) and 3.51 (division) --
% which makes sense: merging two vertices that are already nearly
% coincident barely perturbs the geometry.
param.T_eff_seal = 0.090348;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%%%%%%%%%%%%%%% T1 Intercalation (4-Cell) Parameters %%%%%%%%%%%%%%%%%%%%
% Bulk T1 must be ON to operate in the fluid regime -- Tetley et al. 2019
% identify bulk intercalation as the mechanism that makes fluid tissue
% close wounds faster than solid tissue (their central result).
param.enableT1transitions = true;    % MASTER TOGGLE to turn all T1 swaps on/off
% Wound-margin cells intercalate ~3-4x faster than row-2+ bulk cells in
% Tetley SI Supplementary Fig. 3a-c (row1 vs row2, KS test p<0.0001-0.0079).
% We proxy that row-1 enhancement with a more permissive length tolerance
% at the margin until bulk myosin line-tension fluctuations (their SI
% "Line tension fluctuations", Phase 5 in README_force_based_transitions.md)
% are implemented as the mechanistic driver.
param.enhanceT1atWound = true;     % Toggle for higher intercalation at wound edge
% T1_TOL_bulk = L_T1 / sqrt(A0) from Tetley SI Table 2 (L_T1=0.2um,
% A0=16um^2 -> 0.2/4 = 0.05), rescaled to this model's A0=1.
param.T1_TOL_bulk = 0.05;     % tolerance for bulk tissue
param.T1_TOL_wound = 0.15;    % ~3x bulk tolerance, proxying the row1:row2 ratio above
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%%%%%%%%%%%%%%% Force-Based Transition Parameters %%%%%%%%%%%%%%%%%%%%%%%
param.useForceBasedTransitions = true;
param.k_off0 = 1.0;
% f_beta has no transferable value from the papers (Kaurin-Arroyo's
% f_gamma=sqrt(k*kB*T)~1pN is dimensioned in real pN, not this model's
% dimensionless units) -- their approach is to anchor it at the
% *characteristic* force scale of the system, so we do the same: run the
% model itself, log the tangential/normal force actually seen by edges
% that pass the length pre-condition (Gate 1), and set f_beta to the
% MEDIAN of that population (robust to the heavy right tail from
% near-collapsed edges, unlike the mean). T1 gates on tangential line
% tension, T2 gates on normal pressure -- physically different
% quantities, so they get separate calibrated scales.
%
% Calibration run: 300 steps, fluid regime (p0=4.0), same wound/tissue
% setup as myMainClassical.m, sampling every edge below its T1 length
% gate (tangential) and every wound edge below L_intercalation_thresh
% (normal), each step, regardless of whether it went on to flip.
%   T1 tangential |F|: n=9425, median=22.271433, mean=37.574578 (mean
%     skewed upward by a heavy tail from near-collapsed edges, max~4054;
%     median is the right statistic here, not the mean).
%   T2 normal |pressure|: n=2430, median=0.494935, mean=0.686143.
param.f_beta_T1 = 22.271433;
param.f_beta_T2 = 0.494935;

% Metropolis-style uphill acceptance for the energy gate (T1 in
% checkT1transitions.m, T2 in checkWoundIntercalations.m): instead of a
% strict always-downhill rule, occasionally accept an energy-INCREASING
% swap with probability exp(-dE/T_eff), so the system can escape a
% locally-stuck configuration via fluctuation the way real active/thermal
% systems (and Kaurin-Arroyo's stochastic bond-breaking) can. No
% paper-derivable value -- calibrated the same way as f_beta_T1/f_beta_T2:
% run the model with the strict rule, log dE for every rejected (uphill)
% attempt, anchor T_eff at the median so a "typical" rejected jump has
% roughly e^-1 (~37%) chance of being accepted instead of 0%.
% Calibration run: 600 steps, fluid regime, strict always-downhill rule
% (T_eff unset) so every rejected attempt got logged via
% T1_REJECTED_DE/T2_REJECTED_DE, then took the median of that population.
%   T1 rejected dE: n=994, median=0.038082, mean=0.075429 (right-skewed
%     as usual -- median used, not mean).
param.T_eff_T1 = 0.038082;
% T_eff_T2 recalibrated after redesigning checkWoundIntercalations.m to a
% genuine T1-style neighbor swap (see performWoundMarginT1.m) instead of
% the old merge-and-strip mechanism -- the old calibration (median 0.53)
% no longer applies to the new mechanics.
%   T2 rejected dE: n=7 (small sample -- in this test wound's geometry,
%     most margin vertices are single-real-neighbor "corners" where a
%     4-cell-style swap has no valid partner and is skipped entirely, not
%     rejected -- only a couple of vertices ever reach the energy gate at
%     all), median=1.868404, mean=1.994442.
param.T_eff_T2 = 1.868404;
param.f_beta_length_slope_T1 = 1.0;
param.f_beta_length_slope_T2 = 1.0;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


% %%%%%%%%%%%%%%%%%%%%% NEW: Edge-Specific Purse-String Parameters %%%%%%%%%%%%%%%%%%
param.lambda_purse_string = 4.0;  % Increased for wound closure (was 2.0)
param.t_ramp_purse_string = 50;   % Faster ramp-up (was 100)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%%%%%%%%%%%%%%%%%%% NEW: Wound-Cell Contractility Parameters %%%%%%%%%%%%%%%%%%%%%
% These apply ONLY to wound-edge cells (param.cellIDtoContract).
% Bulk cells are unaffected.
%
%   ka_wound_factor:     multiplier on area stiffness ka for wound cells
%                        > 1 makes wound cells resist area changes more strongly
%                        (drives area contraction toward A0=1)
%
%   contractility_wound: multiplier on perimeter tension (2/rstiff * (P-p0))
%                        > 1 increases the perimeter-tension driving force
%                        (models actomyosin ring contraction at wound margin)
%
param.ka_wound_factor     = 2.0;  % wound cells have 2x area stiffness
param.contractility_wound = 3.0;  % wound cells have 3x perimeter contractility
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%%%%%%%%%%%%%%%%%%% NEW: Wound-Directed Crawling Parameters %%%%%%%%%%%%%%%%%%%%%%
% Models the early, outward-pointing lamellipodial protrusion (OPTL) that
% Trepat et al. 2014 (Nat. Phys.) found dominates the FIRST phase of wound
% closure, before the purse-string ring (IPTL) takes over. Distinct from
% the purse-string: this is an active per-cell force on wound-margin
% cells, directed toward the wound centroid, that decays over time as
% lambda_purse_string ramps up (see t_ramp_purse_string above).
%
%   crawl_force0:   initial magnitude of the crawling force per wound cell
%   tau_decay_crawl: e-folding decay time (in timesteps) of the crawling force
%
param.enableWoundCrawling = true;
param.crawl_force0    = 3.0;
param.tau_decay_crawl = 100;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%%%%%%%%%%%%%%%%%%% Force-Based Closure Strategy (see force_based_closure_strategy.md) %%%%%%%%%%%%%%%%%%%%%%%%%
% Phase F2: wound-hole elastic resistance (Tetley Kw), decays to 0 over
% tau_Kw. Phase F1/F5: purse-string tension, crawl force, and margin
% contractility are all recruited via feedback gated on Kw's decay,
% rather than fixed timers/instant switches. These are first-pass,
% reasoned defaults -- NOT calibrated against logged simulation data the
% way f_beta_T1/f_beta_T2 were; a follow-up Phase-F6 calibration pass
% would set these properly.
param.Kw0    = param.ka * 2;   % initial wound-hole area-elastic modulus, comparable scale to normal cell area stiffness
param.tau_Kw = 2;              % decay time (simulation time units) for Kw -- reasoned to occupy the early fraction of a run, not measured
param.k_recruit = 1.0;         % rate constant for purse-string/contractility recruitment ODEs -- reasoned to match tau_Kw's order of magnitude
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
    % p0* ~ 3.81 is the fluid-solid jamming point (Bi et al. 2015; used by
    % both Tetley SI Fig. S5e-f and README_force_based_transitions.md).
    % Tetley's own shear modulus/viscosity data flattens to ~0 for p0>3.8;
    % 4.0 sits comfortably past that so the tissue is unambiguously fluid.
    param.p0 = 4.0;
    % parm.multFactorForContraction = min(param.multFactorForContraction,1.5); % No longer needed
elseif strcmp(Case,'solid')
    param.p0 = 2;
    % parm.multFactorForContraction = min(param.multFactorForContraction,2.5); % No longer needed
    fprintf(1,'Tissue is now %s\n',Case);
else
    error('Please choose a given tissue fluidity case : "fluid" or "solid".')
end