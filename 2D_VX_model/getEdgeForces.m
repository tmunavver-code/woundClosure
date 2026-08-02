function edgeForceData = getEdgeForces(celldata, param, tstep)
%% getEdgeForces: Compute the force/tension acting along each cell edge
%
% In the classical vertex model, the energy is:
%   E = sum_cells [ ka*(A - A0)^2 + (1/rstiff)*(P - p0)^2 ]
%
% The force contribution along an edge comes from two sources:
%   1. PERIMETER TENSION: The perimeter elasticity term produces a line
%      tension T_peri = 2/rstiff * (P - p0) along each edge. This acts
%      to shorten/lengthen the edge depending on whether P > p0 or P < p0.
%   2. AREA PRESSURE: The area elasticity term produces a pressure-like
%      force normal to each edge. The magnitude is proportional to
%      2*ka*(A - 1).
%
% Additionally, if wound purse-string forces are active, wound-edge
% segments carry an extra contractile tension lambda.
%
% OUTPUT:
%   edgeForceData - struct array (one entry per unique edge) with fields:
%     .v1, .v2         - vertex indices defining the edge
%     .midpoint        - [x, y] midpoint of the edge
%     .edgeVec         - [dx, dy] vector from v1 to v2
%     .length          - scalar length of the edge
%     .tension_peri    - net perimeter tension along the edge (scalar, positive = contractile)
%     .pressure_area   - net area pressure across the edge (scalar, positive = expansive)
%     .tension_purse   - purse-string tension if wound edge, else 0
%     .total_tension   - total effective tension = tension_peri + tension_purse
%     .force_tangent   - [fx, fy] tangential force vector along the edge
%     .force_normal    - [fx, fy] normal (pressure) force vector on the edge
%     .cellIDs         - list of cell IDs sharing this edge
%     .isWoundEdge     - boolean, true if this edge is on the wound margin

Lx = param.Lx;
Ly = param.Ly;
ka = param.ka;
rstiff = param.rstiff;
p0 = param.p0;

%% Step 1: Build a unique edge list with associated cell data
% Use a containers.Map keyed by 'min_max' vertex string to avoid duplicates
edgeMap = containers.Map('KeyType', 'char', 'ValueType', 'any');

for cellID = 1:celldata.nCells
    vertices = celldata.connec{cellID};
    if length(vertices) < 3
        continue;
    end
    
    Acell = celldata.A(cellID);
    Pcell = celldata.P(cellID);
    
    % Perimeter tension for this cell (scalar, same for all edges of this cell)
    T_peri_cell = 2.0 / rstiff * (Pcell - p0);
    
    % Area pressure for this cell
    P_area_cell = 2.0 * ka * (Acell - 1.0);
    
    nVerts = length(vertices);
    for j = 1:nVerts
        v1 = vertices(j);
        if j < nVerts
            v2 = vertices(j + 1);
        else
            v2 = vertices(1);
        end
        
        % Create a canonical key (smaller index first)
        eKey = sprintf('%d_%d', min(v1, v2), max(v1, v2));
        
        if edgeMap.isKey(eKey)
            entry = edgeMap(eKey);
            entry.cellIDs(end + 1) = cellID;
            entry.T_peri_contributions(end + 1) = T_peri_cell;
            entry.P_area_contributions(end + 1) = P_area_cell;
            edgeMap(eKey) = entry;
        else
            entry.v1 = min(v1, v2);
            entry.v2 = max(v1, v2);
            entry.cellIDs = cellID;
            entry.T_peri_contributions = T_peri_cell;
            entry.P_area_contributions = P_area_cell;
            edgeMap(eKey) = entry;
        end
    end
end

%% Step 2: Build wound edge lookup
isWoundEdgeLookup = containers.Map('KeyType', 'char', 'ValueType', 'logical');
if isfield(celldata, 'woundEdges') && ~isempty(celldata.woundEdges)
    for i = 1:size(celldata.woundEdges, 1)
        wv1 = celldata.woundEdges(i, 1);
        wv2 = celldata.woundEdges(i, 2);
        wKey = sprintf('%d_%d', min(wv1, wv2), max(wv1, wv2));
        isWoundEdgeLookup(wKey) = true;
    end
end

%% Step 3: Compute forces for each unique edge
allKeys = edgeMap.keys();
nEdges = length(allKeys);

% Pre-allocate output struct array
edgeForceData(nEdges) = struct('v1', 0, 'v2', 0, 'midpoint', [0 0], ...
    'edgeVec', [0 0], 'length', 0, ...
    'tension_peri', 0, 'pressure_area', 0, 'tension_purse', 0, ...
    'total_tension', 0, ...
    'force_tangent', [0 0], 'force_normal', [0 0], ...
    'cellIDs', [], 'isWoundEdge', false);

% Purse-string tension (with ramp)
if isfield(param, 'lambda_purse_string')
    current_lambda = param.lambda_purse_string * min(1.0, tstep / param.t_ramp_purse_string);
else
    current_lambda = 0;
end

for i = 1:nEdges
    entry = edgeMap(allKeys{i});
    v1 = entry.v1;
    v2 = entry.v2;
    
    % Compute edge vector with periodic boundary handling
    r1 = celldata.r(v1, :);
    r2 = celldata.r(v2, :);
    edgeVec = r2 - r1;
    edgeVec(1) = edgeVec(1) - Lx * round(edgeVec(1) / Lx);
    edgeVec(2) = edgeVec(2) - Ly * round(edgeVec(2) / Ly);
    
    edgeLen = norm(edgeVec);
    midpoint = r1 + 0.5 * edgeVec;
    
    % Wrap midpoint into the box
    midpoint(1) = mod(midpoint(1), Lx);
    midpoint(2) = mod(midpoint(2), Ly);
    
    % Net perimeter tension: sum contributions from all cells sharing this edge
    % Positive = contractile (edge wants to shorten)
    tension_peri = sum(entry.T_peri_contributions);
    
    % Net area pressure: sum contributions from all cells sharing this edge
    % Positive = expansive (cells want to grow)
    pressure_area = sum(entry.P_area_contributions);
    
    % Purse-string tension (only for wound edges)
    eKey = allKeys{i};
    if isWoundEdgeLookup.isKey(eKey)
        tension_purse = current_lambda;
        isWound = true;
    else
        tension_purse = 0;
        isWound = false;
    end
    
    % Total effective line tension along the edge
    total_tension = tension_peri + tension_purse;
    
    % Force vectors
    if edgeLen > 1e-10
        tangent = edgeVec / edgeLen;
        normal = [-tangent(2), tangent(1)]; % 90-degree CCW rotation
        
        % Tangential force (along the edge, contractile if positive)
        force_tangent = total_tension * tangent;
        
        % Normal force (pressure pushing outward from edge midpoint)
        force_normal = pressure_area * edgeLen * normal;
    else
        force_tangent = [0, 0];
        force_normal = [0, 0];
    end
    
    % Store
    edgeForceData(i).v1 = v1;
    edgeForceData(i).v2 = v2;
    edgeForceData(i).midpoint = midpoint;
    edgeForceData(i).edgeVec = edgeVec;
    edgeForceData(i).length = edgeLen;
    edgeForceData(i).tension_peri = tension_peri;
    edgeForceData(i).pressure_area = pressure_area;
    edgeForceData(i).tension_purse = tension_purse;
    edgeForceData(i).total_tension = total_tension;
    edgeForceData(i).force_tangent = force_tangent;
    edgeForceData(i).force_normal = force_normal;
    edgeForceData(i).cellIDs = entry.cellIDs;
    edgeForceData(i).isWoundEdge = isWound;
end

end
