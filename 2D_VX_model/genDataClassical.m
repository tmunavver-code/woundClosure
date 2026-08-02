function celldata = genDataClassical(param)
%getVoronoiTessellation: Generate voronoi tessellation

% For random voronoi
lboundary = [0; param.Lx];
wboundary = [0; param.Ly];
nl = ceil(param.Lx);
nc = ceil(param.Ly);

% Random Sequential addition with contact inhibition
MAXATTEMPT = 10000;
Lx = (lboundary(2) - lboundary(1));
Ly = (wboundary(2) - wboundary(1));
N = nl*nc;

phi_desired = 0.5;	% Desired volume fraction
R_inhibition = sqrt(phi_desired*Lx*Ly/pi/N); % Cell radius

for i = 1:N
    attempt = 0;
    flag = 1;
    while flag == 1 && attempt < MAXATTEMPT
        Xnuclei(i) = lboundary(1) + (lboundary(2) - lboundary(1))*rand();
        Ynuclei(i) = wboundary(1) + (wboundary(2) - wboundary(1))*rand();
        if i > 1
            flag = checkInhibition(i,Xnuclei,Ynuclei,R_inhibition, Lx, Ly);
            attempt = attempt + 1;
        else
            flag = 0;
        end
    end
    
    if attempt == MAXATTEMPT
        fprintf(1,'MAX ATTEMPTS %d reached for particle %d, phi = %2.3f\n',attempt, i, i*pi*R_inhibition^2/Lx/Ly);
        break;
    end
    
    fprintf(1,'Added particle %d in %d attempts, phi = %2.3f\n',i, attempt, i*pi*R_inhibition^2/Lx/Ly);
end

[Xnuclei, Ynuclei] = getTile(Xnuclei,Ynuclei,param.Lx,param.Ly);

r = 1;


%% Obtaining the voronoi tessalation for periodic hex grid
% obtaining the voronoi tessalation
PosNuclei = [Xnuclei Ynuclei];
[V,C] = voronoin(PosNuclei);

% Only take cells within the box (based on cells, not vertices as before)
NotWantedCellIndexX = find((Xnuclei >= lboundary(2))|(Xnuclei < lboundary(1)));
NotWantedCellIndexY = find((Ynuclei > r*wboundary(2))|(Ynuclei < r*wboundary(1)));

nCells = 0;
for i = 1:length(C)
    if find(i == NotWantedCellIndexX)
    elseif find(i == NotWantedCellIndexY)
    else
        nCells = nCells+1; % the number of cells that we keep is the number of closed polygons
        closedPolygonsID(nCells) = i; % stores all cell number that is closed
    end
end

nCells = 0;
for i = 1:length(C)
    if find(i == NotWantedCellIndexX)
    elseif find(i == NotWantedCellIndexY)
    else
        nCells = nCells+1; % the number of cells that we keep is the number of closed polygons
        closedPolygonsID(nCells) = i; % stores all cell number that is closed
    end
    
end

NucleiWithinBox = PosNuclei(closedPolygonsID,:);

% throwing error if Voronoi did not get cells within the desired frame
if nCells == 0
    error('There are no closed cells with the desired domain window.')
end

% Renumber the cells and vertices that are within the calculation domain
[Vsubset, Csubset] = getCellDataWithinDomain(nCells, V, C, closedPolygonsID);

% Get cell barycenter coordinates
cellcenter = zeros(nCells,2);
for cellID = 1:nCells
    cellcenter(cellID,:) = getCellBarycenter(cellID,Csubset,Vsubset);
end

% Rearrange the cell vertices in an anticlockwise order
for cellID = 1:nCells
    Csubset = rearrangeAnticlockwise(cellID,Csubset,Vsubset,cellcenter(cellID,:));
end


% Implement periodicity: removes the image nodes and updates the
% connectivity accordingly (keeping only the nodes inside the box). Periodicity is implemented by calculating the
% forces on one set of nodes, and simply shifting the nodes by one box size
[masterVsubset, masterCsubset] = updateMeshforPeriodicity(Vsubset, Csubset, param);


% assign data to celldata
celldata.nCells  = nCells;
celldata.nMasterVertices = length(masterVsubset);
celldata.r       = masterVsubset;
celldata.connec  = masterCsubset;
celldata.r0      = masterVsubset;             % previous timestep values
celldata.rinit   = masterVsubset;             % previous timestep values
celldata.f       = zeros(size(masterVsubset)); % forces on the vertices

% initialize cell areas and perimeters
celldata.A                      = getCellAreas(celldata,param);
[celldata.P, celldata.EdgeData] = getCellPerimeters(nCells,masterVsubset,masterCsubset,param,0);
celldata.Ainit                  = celldata.A; % Initial cell areas
celldata.verttocell             = getVertexToCellIDMap(celldata.nMasterVertices,nCells,celldata.connec);


%% Obtaing boundary nodes
% Creating a map of the slaves (images) that have been removed
idxA = ismembertol(Vsubset,masterVsubset, 1e-6, 'ByRows', true);
slaveVSubset = Vsubset(idxA == 0, :);
[SlaveToMasterVertexMap, ~] = findMasters(slaveVSubset,masterVsubset,Vsubset,Csubset,celldata,param);

% Master Slave mapping for periodicity
celldata.slave_r          = slaveVSubset;
celldata.StoMmap          = SlaveToMasterVertexMap;

% identifying boundary vertices for stretching experiments
celldata = boundaryVertices(celldata,param);


%% Plotting mesh with cell IDs
figure(1);
[Cnew, Vnew] = getCellDataforPlottingwithoutPeriJumps(nCells, celldata.r, celldata.connec, param);
for cellID = 1:length(closedPolygonsID)
    patch('faces',Cnew{cellID},'vertices',Vnew , ...
        'facecolor','flat','edgecolor',[.2,.2,.2],'CData',cellID)  ;
    barycenter = getCellBarycenter(cellID, Cnew, Vnew);
    text(barycenter(1),barycenter(2),num2str(cellID))
    axis image off;
end
rectangle('Position',[lboundary(1) wboundary(1) (lboundary(2) - lboundary(1)) (wboundary(2) - wboundary(1))])
pause(0.1)

%% (uncomment if needed)

% celldata.slave_r0         = [slaveVSubset];
% celldata.allVertices      = [masterVsubset;slaveVSubset];
% celldata.nAllvertices    = length(celldata.allVertices);
% celldata.augmentedConnec  = SlaveToMasterCellMap;
% [~,celldata.augmentedEdgeData] = getCellPerimeters(nCells,[masterVsubset;slaveVSubset],SlaveToMasterCellMap,param,1);
% celldata.augmentedverttocell             = getVertexToCellIDMap(celldata.nAllvertices,nCells,celldata.augmentedConnec);



% % Plotting network of vertices
% figure(10); hold on;
% rectangle('Position',[lboundary(1) wboundary(1) (lboundary(2) - lboundary(1)) (wboundary(2) - wboundary(1))])
% xlim([-2 10])
% ylim([-2 10])
% voronoi(PosNuclei(:,1),PosNuclei(:,2))
% h = voronoi(NucleiWithinBox(:,1),NucleiWithinBox(:,2));
% set(h,'Color','red', 'LineWidth', 2)
% plot(V(:,1),V(:,2), 'kx', 'MarkerSize', 12)
% plot(Vsubset(:,1),Vsubset(:,2), 'go', 'MarkerSize', 16, 'MarkerFaceColor', 'green')
% plot(masterVsubset(:,1),masterVsubset(:,2), 'bo', 'MarkerSize', 12, 'MarkerFaceColor', 'blue')
% plot(slaveVSubset(:,1),slaveVSubset(:,2), 'bo', 'MarkerSize', 8, 'MarkerFaceColor', 'yellow')


end

