%% Wound set-up only: builds a wounded 30x30 tissue and plots its initial state
addpath(genpath('../common'))
addpath(genpath('./functions2D'))

close all;
clear;
clc;

rng('default');

%% 1 - Initialising the parameters

%%% Choose tissue size
param.Lx       = 30;                 % box length in x
param.Ly       = 30;                 % box length in y

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
param.woundCircular.radius = 1.0;                      % Wound radius

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
param.cellIDtoContract = woundEdgeCells; % The cells that will pull
celldata.woundEdges = woundEdges;


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


%% 3 - Computing
% Storage for a time loop (see myMainClassical.m for the full simulation)

totalsimtime = tic;
final_tstep = param.Nsteps;
energymat = zeros(final_tstep,1);
