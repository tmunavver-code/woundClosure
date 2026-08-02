%% Classical 2D vertex model allowing for T1 transitions
% Simple code to simulate the classical vertex model based on a psuedo-energy
% functional and nodal dissipative forces. The energy minimization also
% allows for topological transitions of type T1.
% Sohan Kale - Adam Ouzeri | September 9 2019
% Modified -  Nikhil Walani October 19 2023

addpath('../common')
addpath('./functions2D')

close all;
clear;
clc;
rng('default');

%% 1 - Initialising the parameters

%%% Choose tissue size
param.Lx       = 15;                 % box length in x
param.Ly       = 15;                 % box length in y

%%% Choose tissue fluidity (comment unwanted)
% Case = 'fluid';
 Case = 'solid';

%%% Choose a scenario (comment unwanted)
 SingleCellContractionParameters;
 %StretchingParameters;
% SelfPropelledParameters;

%% 2 - Generating mesh
celldata = genDataClassical(param);

% plotting initial state
figure(2)
p = plot3Dtissue(celldata.nCells, celldata.r, celldata.connec,param);
rectangle('Position',[0 0 param.Lx param.Ly]);
pause(10)

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
    celldata.f = getVertexForcesClassical(celldata,param,tstep);
    
    %% Update vertex positions while maintaining periodicity
    celldata   = updateVertexPositions(celldata,param);
    
    %% Update cell area and perimeters from new vertex positions
    celldata.A = getCellAreas(celldata,param);
    [celldata.P, celldata.EdgeData] = getCellPerimeters(celldata.nCells,celldata.r,celldata.connec,param,0);
    
    %% Check for T1 transitions and update
    % After performing the T1-swap, edge is relaxed for maxT1relaxsteps
    % steps while keeping all other potential T1-swaps on hold
    [celldata.r,celldata.connec,celldata.EdgeData,celldata.verttocell, T1flagVec,T1relaxstepcountVec,param.nT1]   = checkT1transitions(celldata.nCells,celldata.r,celldata.connec,celldata.EdgeData,celldata.verttocell,param,T1flagVec,T1relaxstepcountVec,0);
    
    
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

