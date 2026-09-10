%% Classical 2D vertex model for wound healing
%
% MODIFIED:
% 1. ADDED: Force visualization plot inside the main time loop.
%
addpath(genpath('../common'))
addpath(genpath('./functions2D'))

close all;
clear;
clc;

rng('default');

%% 1 - Initialising the parameters

%%% Choose tissue size
param.Lx       = 15;                 % box length in x
param.Ly       = 15;                 % box length in y

%%% Choose tissue fluidity
Case = 'fluid'; % 'fluid' or 'solid'

%%% Load base parameters
SingleCellContractionParameters; % Defines ka, rstiff, p0, etc.
% StretchingParameters;
% SelfPropelledParameters;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% --- WOUND CONFIGURATION --- %%%
% Choose the type of wound to create.
% Options: 'circular', 'elliptical', 'none'
param.woundType = 'circular';

% Parameters for a CIRCULAR wound
param.woundCircular.center = [param.Lx/2, param.Ly/2]; % Wound center [x, y]
param.woundCircular.radius = 2.0;                      % Wound radius

% Parameters for an ELLIPTICAL wound
param.woundElliptical.center = [param.Lx/2, param.Ly/2]; % Wound center [x, y]
param.woundElliptical.x_axis = 8.0;                      % Semi-axis in x-direction
param.woundElliptical.y_axis = 1.0;                      % Semi-axis in y-direction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 2 - Generating mesh and Creating Wound
celldata = genDataClassical(param);


%%% Create the specified wound by removing cells
switch param.woundType
    case 'circular'
        [celldata, ~] = createWound(celldata, 'circular', param.woundCircular, param);
    case 'elliptical'
        [celldata, ~] = createWound(celldata, 'elliptical', param.woundElliptical, param);
    case 'none'
        disp('No wound will be created.');
    otherwise
        disp('Invalid wound type specified. No wound created.');
end


%%% Identify wound edge cells for purse-string contraction
[woundEdgeCells, woundEdges] = findWoundEdgeCells(celldata);
param.cellIDtoContract = woundEdgeCells; % <-- THIS IS STILL NEEDED for plot3Dtissue coloring
celldata.woundEdges = woundEdges; % <-- THIS IS NOW THE IMPORTANT PART
celldata.cellIDtoContract = woundEdgeCells;

%%% Force-based closure state (see force_based_closure_strategy.md)
% Kw0/tau_Kw/k_recruit come from SingleCellContractionParameters.m.
[woundLoopVerts0, woundLoopValid0] = getOrderedWoundLoop(celldata.woundEdges);
if woundLoopValid0
    celldata_tmp0 = celldata;
    celldata_tmp0.connec{celldata.nCells+1} = woundLoopVerts0;
    celldata.woundArea0 = getPolygonalCellArea(celldata_tmp0, celldata.nCells+1, param);
else
    celldata.woundArea0 = 0; % no valid single wound loop -- Kw force stays off (see getOrderedWoundLoop.m)
end
celldata.Kw_current = param.Kw0;
celldata.lambda_current = 0;
celldata.crawl_current = 0;
celldata.ka_wound_factor_current = 1;
celldata.contractility_wound_current = 1;


% Plotting initial state with the wound
figure(2)
plot3Dtissue(celldata.nCells, celldata.r, celldata.connec,param);
hold on;
% Optional: Highlight the purse-string cells
p_wound = plot3Dtissue(length(woundEdgeCells), celldata.r, celldata.connec(woundEdgeCells), param);
set(p_wound(isgraphics(p_wound)), 'FaceColor', 'magenta', 'FaceAlpha', 0.5);
title('Initial State with Wound');
rectangle('Position',[0 0 param.Lx param.Ly]);
hold off;
pause(5)


%% 3 - Computing
totalsimtime = tic;

energymat = zeros(param.Nsteps,1);
timemat   = zeros(param.Nsteps,1);
T1flagVec    = zeros(celldata.nMasterVertices,1);
T1relaxstepcountVec = zeros(size(T1flagVec));
T1ForceLog = [];

Coordinates = zeros(celldata.nMasterVertices,2*param.Nsteps);
Connectivity = cell(celldata.nCells,param.Nsteps);

%% Time loop
for tstep = 1:param.Nsteps
    
    %% Applying stretch (does nothing when not wanted)
    applyStretchX;
    
    
    %% initial energy
    [energymat0] = getTissueEnergyClassical(celldata,param);

    %% Force-based closure state updates (see force_based_closure_strategy.md)
    % Kw decays analytically; purse-string tension, crawl force, and
    % margin contractility are all recruited via feedback ODEs gated on
    % how much of that decay has happened, instead of independent fixed
    % timers.
    if isfield(celldata, 'woundArea0') && celldata.woundArea0 > 0
        celldata.Kw_current = param.Kw0 * exp(-tstep * param.deltat / param.tau_Kw);
        kw_fraction_decayed = 1 - celldata.Kw_current / param.Kw0;
    else
        kw_fraction_decayed = 1; % no wound (or invalid loop) -- nothing to gate on
    end

    target_lambda = param.lambda_purse_string * kw_fraction_decayed;
    celldata.lambda_current = celldata.lambda_current + param.deltat * param.k_recruit * (target_lambda - celldata.lambda_current);

    if param.lambda_purse_string > 0
        celldata.crawl_current = param.crawl_force0 * max(0, 1 - celldata.lambda_current / param.lambda_purse_string);
    else
        celldata.crawl_current = param.crawl_force0;
    end

    target_ka_factor = 1 + (param.ka_wound_factor - 1) * kw_fraction_decayed;
    celldata.ka_wound_factor_current = celldata.ka_wound_factor_current + param.deltat * param.k_recruit * (target_ka_factor - celldata.ka_wound_factor_current);

    target_contractility_factor = 1 + (param.contractility_wound - 1) * kw_fraction_decayed;
    celldata.contractility_wound_current = celldata.contractility_wound_current + param.deltat * param.k_recruit * (target_contractility_factor - celldata.contractility_wound_current);

    %% Get forces based on the current configuration
    % This function now includes elastic forces + edge-specific purse-string
    celldata.f = getVertexForcesClassical(celldata,param,tstep);


    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ---- VISUAL REPRESENTATION OF FORCES ---
    % Moved here (was after the position update) so the plotted force
    % vectors are paired with the positions that actually produced them.
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if param.plotForces && mod(tstep, param.plotForces_Interval) == 0

        figure(param.plotForces_FigHandle); % Plot in a dedicated window
        clf;

        % Plot the tissue first
        plot3Dtissue(celldata.nCells, celldata.r, celldata.connec, param);
        hold on;

        % Highlight the wound-edge cells
        if ~isempty(woundEdgeCells)
            p_wound = plot3Dtissue(length(woundEdgeCells), celldata.r, celldata.connec(woundEdgeCells), param);
            set(p_wound(isgraphics(p_wound)), 'FaceColor', 'magenta', 'FaceAlpha', 0.5);
        end

        % Get vertex and force data
        r_x = celldata.r(:, 1);
        r_y = celldata.r(:, 2);
        f_x = celldata.f(:, 1);
        f_y = celldata.f(:, 2);

        % Plot force vectors using quiver
        quiver(r_x, r_y, f_x, f_y, 'r', 'AutoScaleFactor', param.plotForces_Scale);

        title(['Forces at Timestep: ', num2str(tstep)]);
        rectangle('Position',[0 0 param.Lx param.Ly]);
        axis equal; % Ensure aspect ratio is correct
        hold off;

        drawnow; % Force MATLAB to render the plot
    end
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %% Update vertex positions while maintaining periodicity
    celldata   = updateVertexPositions(celldata,param);

    %% Update cell area and perimeters from new vertex positions
    celldata.A = getCellAreas(celldata,param);
    [celldata.P, celldata.EdgeData] = getCellPerimeters(celldata.nCells,celldata.r,celldata.connec,param,0);


    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute per-edge forces ONCE on this step's fresh, self-consistent
    % geometry and reuse for both transition checks below -- avoids gating
    % on forces computed before this step's position update.
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    edgeForceData = getEdgeForces(celldata, param, tstep);


    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % --- 3-CELL WOUND INTERCALATION (PRIORITY) ---
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if param.enableWoundIntercalations

        [celldata, T1flagVec, intercalation_happened, param] = checkWoundIntercalations(celldata, param, T1flagVec, edgeForceData, tstep);

        if intercalation_happened
            celldata.A = getCellAreas(celldata,param);
            [celldata.P, celldata.EdgeData] = getCellPerimeters(celldata.nCells,celldata.r,celldata.connec,param,0);
            % Geometry changed -- edgeForceData is now stale for the
            % touched edges. Recompute before the T1 check below.
            edgeForceData = getEdgeForces(celldata, param, tstep);
        end
    end
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %% Check for T1 transitions and update
    % This handles the 4-cell "fluidity" swaps. Single call per step now
    % (the previous duplicate pre/post-update call used a stale geometry
    % snapshot for the first one and has been removed).
    nT1_before = param.nT1;
    if param.enableT1transitions
        [celldata.r,celldata.connec,celldata.EdgeData,celldata.verttocell, T1flagVec,T1relaxstepcountVec,param.nT1, T1ForceLog]   = ...
            checkT1transitions(celldata.nCells,celldata.r,celldata.connec,celldata.EdgeData,celldata.verttocell,param,T1flagVec,T1relaxstepcountVec,0, edgeForceData, tstep, T1ForceLog, celldata);
    end
    if param.nT1 > nT1_before
        % Bulk T1 rewired connectivity -- refresh the wound margin (see
        % refreshWoundMargin.m).
        [celldata, param] = refreshWoundMargin(celldata, param);
    end

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % --- CELL DIVISION (force-based, third topological event type) ---
    % See checkCellDivision.m / performCellDivision.m / force_based_closure_
    % strategy.md. Adds new vertices (nMasterVertices grows), so every
    % other per-vertex array sized at the OLD nMasterVertices must grow
    % with it, or later indexing into celldata.r/celldata.f crashes.
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nMasterVerticesBefore = celldata.nMasterVertices;
    nCellsBefore = celldata.nCells;
    [celldata, param, division_happened] = checkCellDivision(celldata, param);
    if division_happened
        nNewVerts = celldata.nMasterVertices - nMasterVerticesBefore;
        if nNewVerts > 0
            T1flagVec(end+1:end+nNewVerts) = 0;
            T1relaxstepcountVec(end+1:end+nNewVerts) = 0;
            celldata.f(end+1:end+nNewVerts, :) = 0;
            celldata.r0(end+1:end+nNewVerts, :) = celldata.r(end-nNewVerts+1:end, :);
            Coordinates(end+1:end+nNewVerts, :) = 0;
        end
        nNewCells = celldata.nCells - nCellsBefore;
        if nNewCells > 0
            Connectivity(end+1:end+nNewCells, :) = {[]};
        end
        celldata.A = getCellAreas(celldata,param);
        [celldata.P, celldata.EdgeData] = getCellPerimeters(celldata.nCells,celldata.r,celldata.connec,param,0);
        % Division rewired connectivity and added cells -- refresh margin.
        [celldata, param] = refreshWoundMargin(celldata, param);
    end

    %% Wound sealing (force + energy gated) -- retires cells from the
    % margin as opposing sides come into contact. See checkWoundSealing.m.
    [celldata, param, seal_happened] = checkWoundSealing(celldata, param, tstep);

    %% Storing current state for plotting purposes
    Coordinates(:,2*tstep - 1:2*tstep) = celldata.r;
    Connectivity(:,tstep)= celldata.connec;
    
    %% Update time if no more T1 possible
    param.Tsim = param.Tsim + param.deltat;
    celldata.r0 = celldata.r;
    
    %% Update status
    energymat(tstep) = getTissueEnergyClassical(celldata,param);
    timemat(tstep)   = param.Tsim;
    if tstep == 1
        relenergychange  = (energymat(tstep) - energymat0)/ energymat0;
    else
        relenergychange  = (energymat(tstep) - energymat(tstep-1))/ energymat(tstep-1);
    end
    
    PrintInformation;
end

%% Post-Processing
figure(3)
PlotTissueEvolution;

% Plotting evolution of the energy
figure(4); hold on; plot(timemat, energymat, 'LineWidth', 3.0); title('Energy vs time')

%% Save Results
if ~exist('results', 'dir')
    mkdir('results');
end
save('results/myMainClassical_workspace.mat');
if ishandle(2), saveas(figure(2), 'results/myMainClassical_initial_state.png'); end
if ishandle(3), saveas(figure(3), 'results/myMainClassical_tissue_evolution.png'); end
if ishandle(4), saveas(figure(4), 'results/myMainClassical_energy_vs_time.png'); end
if ishandle(10), saveas(figure(10), 'results/myMainClassical_forces.png'); end