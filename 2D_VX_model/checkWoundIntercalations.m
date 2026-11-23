function [celldata, T1flagVec, intercalation_happened, param, eventInfo] = ...
    checkWoundIntercalations(celldata, param, T1flagVec, edgeForceData, tstep)
% CHECKWOUNDINTERCALATIONS Performs validated wound-margin neighbor exchanges.
% True three-cell handoffs are prioritized over corner simplification.

intercalation_happened = false;
eventInfo = emptyEventInfo(tstep);

if ~isfield(celldata, 'woundEdges') || isempty(celldata.woundEdges)
    return;
end

coordinates = celldata.r;
connectivity = celldata.connec;
verttocell = celldata.verttocell;
woundEdges = celldata.woundEdges;
woundCells = param.cellIDtoContract;
oldBoundaryCount = numel(unique(woundCells));
oldDegrees = cellfun(@(x) numel(unique(x)), verttocell);
oldMaxDegree = max(oldDegrees);

% Target boundary population shrinks with the measured wound perimeter
% (lower bound param.minWoundBoundaryCells).
currentPerimeter = woundPerimeter(coordinates, woundEdges, param);
if isfield(celldata, 'woundPerimeter0') && celldata.woundPerimeter0 > 0
    perimeterRatio = min(max(currentPerimeter / celldata.woundPerimeter0, 0), 1);
else
    perimeterRatio = 1;
end
if isfield(param, 'initialWoundBoundaryCells')
    initialBoundaryCount = param.initialWoundBoundaryCells;
else
    initialBoundaryCount = oldBoundaryCount;
end
if isfield(param, 'minWoundBoundaryCells')
    minimumBoundaryCount = param.minWoundBoundaryCells;
else
    minimumBoundaryCount = 5;
end
targetBoundaryCount = max(minimumBoundaryCount, ...
    ceil(initialBoundaryCount * perimeterRatio));

if isfield(param, 'bellProbabilityThresholdT2')
    probabilityThresholdMax = param.bellProbabilityThresholdT2;
else
    probabilityThresholdMax = param.bellProbabilityThreshold;
end
if isfield(param, 'bellProbabilityThresholdT2Min')
    probabilityThresholdMin = param.bellProbabilityThresholdT2Min;
else
    probabilityThresholdMin = probabilityThresholdMax;
end
probabilityThreshold = probabilityThresholdMin + ...
    (probabilityThresholdMax - probabilityThresholdMin) * perimeterRatio;

% A handoff is allowed only while the margin holds more cells than the target.
allowHandoff = oldBoundaryCount > targetBoundaryCount;
if allowHandoff
    % Crowded margin: use the lower Bell threshold.
    probabilityThreshold = probabilityThresholdMin;
end
eventInfo.probabilityThreshold = probabilityThreshold;
eventInfo.perimeterRatio = perimeterRatio;
eventInfo.targetBoundaryCount = targetBoundaryCount;

candidates = struct('edgeIndex', {}, 'v1', {}, 'v2', {}, 'cellA', {}, ...
    'cellB', {}, 'cellC', {}, 'mode', {}, 'modeRank', {}, 'length', {}, ...
    'lengthRatio', {}, 'forceDrive', {}, 'probability', {}, 'hasForce', {});

% Classify every boundary edge before choosing an event.
for edgeIndex = 1:size(woundEdges, 1)
    v1 = woundEdges(edgeIndex, 1);
    v2 = woundEdges(edgeIndex, 2);
    if T1flagVec(v1) == 1 || T1flagVec(v2) == 1
        eventInfo.rejectedRelaxation = eventInfo.rejectedRelaxation + 1;
        continue;
    end

    cellsAtV1 = verttocell{v1};
    cellsAtV2 = verttocell{v2};
    cellAList = intersect(intersect(cellsAtV1, cellsAtV2), woundCells);
    if numel(cellAList) ~= 1
        eventInfo.rejectedTopology = eventInfo.rejectedTopology + 1;
        continue;
    end
    cellA = cellAList(1);
    if numel(connectivity{cellA}) < 4
        eventInfo.rejectedCellGeometry = eventInfo.rejectedCellGeometry + 1;
        continue;
    end

    cellsB = setdiff(cellsAtV1, cellA);
    cellsC = setdiff(cellsAtV2, cellA);
    cellB = NaN;
    cellC = NaN;
    if numel(cellsB) == 1 && numel(cellsC) == 1 && cellsB(1) ~= cellsC(1)
        mode = 'handoff';
        modeRank = 1;
        cellB = cellsB(1);
        cellC = cellsC(1);
    elseif isempty(cellsB) && numel(cellsC) == 1
        mode = 'trim_v1';
        modeRank = 2;
        cellC = cellsC(1);
    elseif isempty(cellsC) && numel(cellsB) == 1
        mode = 'trim_v2';
        modeRank = 2;
        cellB = cellsB(1);
    else
        eventInfo.rejectedTopology = eventInfo.rejectedTopology + 1;
        continue;
    end

    d = coordinates(v2,:) - coordinates(v1,:);
    d(1) = d(1) - param.Lx * round(d(1) / param.Lx);
    d(2) = d(2) - param.Ly * round(d(2) / param.Ly);
    edgeLength = norm(d);
    localEdgeScale = celldata.P(cellA) / numel(connectivity{cellA});
    lengthRatio = edgeLength / max(localEdgeScale, 1e-12);

    forceDrive = NaN;
    probability = 1;
    hasForce = false;
    if nargin >= 4 && ~isempty(edgeForceData)
        vMin = min(v1, v2);
        vMax = max(v1, v2);
        forceIndex = find([edgeForceData.v1] == vMin & ...
            [edgeForceData.v2] == vMax, 1);
        if ~isempty(forceIndex)
            hasForce = true;
            % Compressive area pressure drives loss of the wound contact.
            forceDrive = max(-edgeForceData(forceIndex).pressure_area, 0);
            fBeta = max(param.f_beta_T2 - param.f_beta_length_slope_T2 * ...
                max(0, param.L_intercalation_thresh - edgeLength), 0.01);
            kOff = param.k_off0 * exp(forceDrive / fBeta);
            probability = 1 - exp(-param.deltat * kOff);
        end
    end

    candidate.edgeIndex = edgeIndex;
    candidate.v1 = v1;
    candidate.v2 = v2;
    candidate.cellA = cellA;
    candidate.cellB = cellB;
    candidate.cellC = cellC;
    candidate.mode = mode;
    candidate.modeRank = modeRank;
    candidate.length = edgeLength;
    candidate.lengthRatio = lengthRatio;
    candidate.forceDrive = forceDrive;
    candidate.probability = probability;
    candidate.hasForce = hasForce;
    candidates(end+1) = candidate; %#ok<AGROW>
end

eventInfo.candidateCount = numel(candidates);
if isempty(candidates)
    return;
end

sortMatrix = [[candidates.modeRank]', [candidates.lengthRatio]'];
[~, order] = sortrows(sortMatrix, [1 2]);
E_old = getTissueEnergyClassical(celldata, param);

for orderIndex = 1:numel(order)
    candidate = candidates(order(orderIndex));
    eventInfo.consideredCount = eventInfo.consideredCount + 1;

    if strcmp(candidate.mode, 'handoff') && ~allowHandoff
        eventInfo.rejectedBoundaryTarget = ...
            eventInfo.rejectedBoundaryTarget + 1;
        continue;
    end

    % Trigger: a short wound contact, or a crowded margin (shortest handoff).
    perimeterCrowdingTrigger = strcmp(candidate.mode, 'handoff') && ...
        allowHandoff;
    lengthEligible = candidate.length < param.L_intercalation_thresh || ...
        candidate.lengthRatio < param.L_intercalation_ratio_thresh || ...
        perimeterCrowdingTrigger;
    if ~lengthEligible
        eventInfo.rejectedLength = eventInfo.rejectedLength + 1;
        continue;
    end
    if isfield(param, 'useForceBasedTransitions') && ...
            param.useForceBasedTransitions
        if ~candidate.hasForce
            eventInfo.rejectedMissingForce = eventInfo.rejectedMissingForce + 1;
            continue;
        elseif candidate.probability < probabilityThreshold
            eventInfo.rejectedProbability = eventInfo.rejectedProbability + 1;
            continue;
        end
    end

    coordsTrial = coordinates;
    connecTrial = connectivity;
    switch candidate.mode
        case 'handoff'
            [coordsTrial, connecTrial, ~] = performWoundMarginT1( ...
                candidate.cellA, candidate.cellB, candidate.cellC, ...
                coordsTrial, connecTrial, verttocell, candidate.v1, ...
                candidate.v2, param);
        case 'trim_v1'
            connecTrial = trimWoundCorner(candidate.cellA, connecTrial, candidate.v1);
        case 'trim_v2'
            connecTrial = trimWoundCorner(candidate.cellA, connecTrial, candidate.v2);
    end

    trial = celldata;
    trial.r = coordsTrial;
    trial.connec = connecTrial;
    trial.verttocell = getVertexToCellIDMap(size(coordsTrial,1), ...
        trial.nCells, connecTrial);
    trial.A = getCellAreas(trial, param);
    [trial.P, trial.EdgeData] = getCellPerimeters(trial.nCells, trial.r, ...
        trial.connec, param, 0);

    [newWoundCells, newWoundEdges] = findWoundEdgeCells(trial);
    trial.woundEdges = newWoundEdges;
    trial.cellIDtoContract = newWoundCells;
    trialParam = param;
    trialParam.cellIDtoContract = newWoundCells;

    [topologyValid, ~, trialMaxDegree] = validateTissueTopology(trial);
    [~, loopValid] = getOrderedWoundLoop(newWoundEdges);
    if isempty(newWoundEdges)
        loopValid = true;
    end
    allowedMaxDegree = max(oldMaxDegree, param.maxSealVertexDegree);
    positiveAreas = trial.A(trial.A > 0);
    if ~topologyValid || ~loopValid || trialMaxDegree > allowedMaxDegree || ...
            any(~isfinite(trial.A)) || isempty(positiveAreas) || ...
            min(positiveAreas) <= 1e-10
        eventInfo.rejectedValidation = eventInfo.rejectedValidation + 1;
        continue;
    end

    if strcmp(candidate.mode, 'handoff') && ...
            (ismember(candidate.cellA, newWoundCells) || ...
            numel(newWoundCells) >= oldBoundaryCount)
        eventInfo.rejectedNoHandoff = eventInfo.rejectedNoHandoff + 1;
        continue;
    end

    E_new = getTissueEnergyClassical(trial, trialParam);
    dE = E_new - E_old;
    if dE < 0
        accept = true;
    elseif isfield(param, 'T_eff_T2') && param.T_eff_T2 > 0
        accept = rand() < exp(-dE / param.T_eff_T2);
    else
        accept = false;
    end
    if ~accept
        eventInfo.rejectedEnergy = eventInfo.rejectedEnergy + 1;
        continue;
    end

    celldata = trial;
    param = trialParam;
    T1flagVec(candidate.v1) = 1;
    T1flagVec(candidate.v2) = 1;
    intercalation_happened = true;
    eventInfo.happened = true;
    eventInfo.mode = candidate.mode;
    eventInfo.cellID = candidate.cellA;
    eventInfo.v1 = candidate.v1;
    eventInfo.v2 = candidate.v2;
    eventInfo.edgeLength = candidate.length;
    eventInfo.lengthRatio = candidate.lengthRatio;
    eventInfo.forceDrive = candidate.forceDrive;
    eventInfo.bellProbability = candidate.probability;
    eventInfo.energyChange = dE;
    eventInfo.boundaryCellsBefore = oldBoundaryCount;
    eventInfo.boundaryCellsAfter = numel(newWoundCells);

    if strcmp(candidate.mode, 'handoff')
        param.nWoundIntercalations = param.nWoundIntercalations + 1;
        fprintf('WOUND_HANDOFF step=%d cell=%d boundary=%d->%d P=%.5f dE=%.5f\n', ...
            tstep, candidate.cellA, oldBoundaryCount, numel(newWoundCells), ...
            candidate.probability, dE);
    else
        param.nWoundCornerTrims = param.nWoundCornerTrims + 1;
        fprintf('WOUND_TRIM step=%d cell=%d mode=%s P=%.5f dE=%.5f\n', ...
            tstep, candidate.cellA, candidate.mode, candidate.probability, dE);
    end
    return;
end

end

function eventInfo = emptyEventInfo(tstep)
eventInfo = struct('timestep', tstep, 'happened', false, 'mode', '', ...
    'cellID', NaN, 'v1', NaN, 'v2', NaN, 'edgeLength', NaN, ...
    'lengthRatio', NaN, 'forceDrive', NaN, 'bellProbability', NaN, ...
    'probabilityThreshold', NaN, 'energyChange', NaN, ...
    'perimeterRatio', NaN, 'targetBoundaryCount', NaN, ...
    'boundaryCellsBefore', NaN, 'boundaryCellsAfter', NaN, ...
    'candidateCount', 0, 'consideredCount', 0, 'rejectedRelaxation', 0, ...
    'rejectedTopology', 0, 'rejectedCellGeometry', 0, 'rejectedLength', 0, ...
    'rejectedMissingForce', 0, 'rejectedProbability', 0, ...
    'rejectedBoundaryTarget', 0, 'rejectedValidation', 0, ...
    'rejectedNoHandoff', 0, 'rejectedEnergy', 0);
end

function perimeter = woundPerimeter(coordinates, woundEdges, param)
perimeter = 0;
for edgeIndex = 1:size(woundEdges, 1)
    d = coordinates(woundEdges(edgeIndex,2),:) - ...
        coordinates(woundEdges(edgeIndex,1),:);
    d(1) = d(1) - param.Lx * round(d(1) / param.Lx);
    d(2) = d(2) - param.Ly * round(d(2) / param.Ly);
    perimeter = perimeter + norm(d);
end
end
