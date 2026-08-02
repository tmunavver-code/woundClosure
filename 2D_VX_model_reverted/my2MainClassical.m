%% Classical 2D vertex model for wound healing
%
% MODIFIED:
% 1. ADDED: Force visualization plot inside the main time loop.
% 2. ADDED: Call to 3-cell intercalation function (checkWoundIntercalations)
% 3. ADDED: Logic to save all plots to a unique, timestamped folder
%    with detailed run parameters in the name.
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
param.woundCircular.radius = 1.5;                      % Wound radius

% Parameters for an ELLIPTICAL wound
param.woundElliptical.center = [param.Lx/2, param.Ly/2]; % Wound center [x, y]
param.woundElliptical.x_axis = 2.0;                      % Semi-axis in x-direction
param.woundElliptical.y_axis = 1.0;                      % Semi-axis in y-direction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- NEW: CREATE SAVE DIRECTORY (with detailed name) ---
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- 1. Build Wound String ---
if strcmp(param.woundType, 'circular')
    woundStr = sprintf('Circ_r%.1f', param.woundCircular.radius);
elseif strcmp(param.woundType, 'elliptical')
    woundStr = sprintf('Ellip_%.1fx%.1f', param.woundElliptical.x_axis, param.woundElliptical.y_axis);
else
    woundStr = 'NoWound';
end

% --- 2. Build Intercalation String ---
if param.enableWoundIntercalations
    intercalStr = 'Intercal_On';
else
    intercalStr = 'Intercal_Off';
end

% --- 3. Build T1 String ---
if param.enableT1transitions
    t1Str = 'T1_On';
else
    t1Str = 'T1_Off';
end

% --- 4. Build Contractility String ---
lambdaStr = sprintf('L_%.2f', param.lambda_purse_string);

% --- 5. Build Full Directory Name ---
dateStr = datestr(now, 'yyyy-mm-dd_HHMMSS');
saveDirName = sprintf('%s_%s_%s_%s_%s_%s', ...
    dateStr, ...
    Case, ...
    woundStr, ...
    intercalStr, ...
    t1Str, ...
    lambdaStr);

% --- 6. Create Directory ---
mkdir(saveDirName);
fprintf('Saving plots to folder: %s\n', saveDirName);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
celldata.cellIDtoContract = woundEdgeCells; % Store in celldata as well


% Plotting initial state with the wound
figure(2)
plot3Dtissue(celldata.nCells, celldata.r, celldata.connec,param);
hold on;
% Optional: Highlight the purse-string cells
p_wound = plot3Dtissue(length(woundEdgeCells), celldata.r, celldata.connec(woundEdgeCells), param);
set(p_wound, 'FaceColor', 'magenta', 'FaceAlpha', 0.5);
title('Initial State with Wound');
rectangle('Position',[0 0 param.Lx param.Ly]);
hold off;
pause(5)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- NEW: SAVE FIGURE 2 ---
fig2_path = fullfile(saveDirName, '01_Initial_State.png');
saveas(figure(2), fig2_path);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 3 - Computing
totalsimtime = tic;

energymat = zeros(param.Nsteps,1);
timemat   = zeros(param.Nsteps,1);
T1flagVec    = zeros(celldata.nMasterVertices,1);
T1relaxstepcountVec = zeros(size(T1flagVec));

Coordinates = zeros(celldata.nMasterVertices,2*param.Nsteps);
Connectivity = cell(celldata.nCells,param.Nsteps);

%% Time loop
for tstep = 1:param.Nsteps
    
    %% Applying stretch (does nothing when not wanted)
    applyStretchX;
    
    
    %% initial energy
    [energymat0] = getTissueEnergyClassical(celldata,param);
    
    %% Get forces based on the current configuration
    % This function now includes elastic forces + edge-specific purse-string
    celldata.f = getVertexForcesClassical(celldata,param,tstep);
    
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % --- NEW: PRIORITY 2 - VISUAL REPRESENTATION OF FORCES ---
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if param.plotForces && mod(tstep, param.plotForces_Interval) == 0
        
        figure(param.plotForces_FigHandle); % Plot in a dedicated window
        clf; 
        
        % Plot the tissue first
        plot3Dtissue(celldata.nCells, celldata.r, celldata.connec, param);
        hold on;
        
        % Highlight the wound-edge cells
        if ~isempty(param.cellIDtoContract)
            p_wound = plot3Dtissue(length(param.cellIDtoContract), celldata.r, celldata.connec(param.cellIDtoContract), param);
            set(p_wound, 'FaceColor', 'magenta', 'FaceAlpha', 0.5);
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
        
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % --- NEW: SAVE FIGURE 10 (Force Plot) ---
        fig10_filename = sprintf('02_Forces_Tstep_%05d.png', tstep); % %05d pads with zeros
        fig10_path = fullfile(saveDirName, fig10_filename);
        saveas(figure(param.plotForces_FigHandle), fig10_path);
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % --- END OF NEW SECTION ---
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %% Update vertex positions while maintaining periodicity
    celldata   = updateVertexPositions(celldata,param);
    
    %% Update cell area and perimeters from new vertex positions
    celldata.A = getCellAreas(celldata,param);
    [celldata.P, celldata.EdgeData] = getCellPerimeters(celldata.nCells,celldata.r,celldata.connec,param,0);
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % --- NEW: 3-CELL WOUND INTERCALATION (PRIORITY) ---
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This event changes the topology, so it must run first.
    if param.enableWoundIntercalations
        
        [celldata, T1flagVec, intercalation_happened] = checkWoundIntercalations(celldata, param, T1flagVec);
        
        % If an event occurred, the geometry is invalid. We MUST recalculate
        % Areas and Perimeters before proceeding to T1 checks or force calcs.
        if intercalation_happened
            celldata.A = getCellAreas(celldata,param);
            [celldata.P, celldata.EdgeData] = getCellPerimeters(celldata.nCells,celldata.r,celldata.connec,param,0);
        end
    end
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    %% Check for T1 transitions and update
    % This handles the 4-cell "fluidity" swaps
    if param.enableT1transitions
        [celldata.r,celldata.connec,celldata.EdgeData,celldata.verttocell, T1flagVec,T1relaxstepcountVec,param.nT1]   = checkT1transitions(celldata.nCells,celldata.r,celldata.connec,celldata.EdgeData,celldata.verttocell,param,T1flagVec,T1relaxstepcountVec,0);
    end
    
    
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
fprintf('Generating and saving final plots (hidden)...\n');

% --- Figure 3: Tissue Evolution (Hidden) ---
fig3 = figure('visible', 'off'); % Create a hidden figure
PlotTissueEvolution; % This script will plot to the active (hidden) figure
fig3_path = fullfile(saveDirName, '03_Tissue_Evolution.png');
saveas(fig3, fig3_path);
close(fig3); % Close the hidden figure

% --- Figure 4: Energy vs Time (Hidden) ---
fig4 = figure('visible', 'off'); % Create another hidden figure
plot(timemat, energymat, 'LineWidth', 3.0);
hold on;
title('Energy vs time');
fig4_path = fullfile(saveDirName, '04_Energy_vs_Time.png');
saveas(fig4, fig4_path);
close(fig4); % Close the hidden figure

fprintf('All plots saved to folder: %s\n', saveDirName);