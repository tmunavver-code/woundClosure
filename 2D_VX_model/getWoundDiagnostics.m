function diagnostics = getWoundDiagnostics(celldata, param)
% GETWOUNDDIAGNOSTICS Summarizes wound geometry and topology for one step.

[boundaryCells, woundEdges] = findWoundEdgeCells(celldata);
diagnostics.boundaryCellCount = numel(boundaryCells);
diagnostics.boundaryEdgeCount = size(woundEdges, 1);
diagnostics.boundaryVertexCount = numel(unique(woundEdges(:)));
diagnostics.areaDeficitRaw = param.Lx * param.Ly - sum(celldata.A);
diagnostics.areaConservation = max(0, diagnostics.areaDeficitRaw);

diagnostics.woundPerimeter = 0;
for edgeID = 1:size(woundEdges, 1)
    d = celldata.r(woundEdges(edgeID,2),:) - celldata.r(woundEdges(edgeID,1),:);
    d(1) = d(1) - param.Lx * round(d(1) / param.Lx);
    d(2) = d(2) - param.Ly * round(d(2) / param.Ly);
    diagnostics.woundPerimeter = diagnostics.woundPerimeter + norm(d);
end

[loopVertices, loopValid] = getOrderedWoundLoop(woundEdges);
diagnostics.loopValid = loopValid || isempty(woundEdges);
diagnostics.loopArea = NaN;
if isempty(woundEdges)
    diagnostics.loopArea = 0;
elseif loopValid
    coords = zeros(numel(loopVertices), 2);
    coords(1,:) = celldata.r(loopVertices(1),:);
    for k = 2:numel(loopVertices)
        d = celldata.r(loopVertices(k),:) - celldata.r(loopVertices(k-1),:);
        d(1) = d(1) - param.Lx * round(d(1) / param.Lx);
        d(2) = d(2) - param.Ly * round(d(2) / param.Ly);
        coords(k,:) = coords(k-1,:) + d;
    end
    x = coords(:,1);
    y = coords(:,2);
    diagnostics.loopArea = 0.5 * abs(sum(x .* circshift(y,-1) - y .* circshift(x,-1)));
end

degrees = cellfun(@(x) numel(unique(x)), celldata.verttocell);
diagnostics.maxVertexDegree = max(degrees);
if isempty(woundEdges)
    diagnostics.maxWoundVertexDegree = 0;
else
    diagnostics.maxWoundVertexDegree = max(degrees(unique(woundEdges(:))));
end

[diagnostics.topologyValid, diagnostics.topologyReason, ~] = ...
    validateTissueTopology(celldata);
positiveAreas = celldata.A(celldata.A > 0);
if isempty(positiveAreas)
    diagnostics.minPositiveCellArea = NaN;
else
    diagnostics.minPositiveCellArea = min(positiveAreas);
end

diagnostics.minEdgeLength = Inf;
for cellID = 1:celldata.nCells
    vertices = celldata.connec{cellID};
    for j = 1:numel(vertices)
        nextJ = mod(j, numel(vertices)) + 1;
        d = celldata.r(vertices(nextJ),:) - celldata.r(vertices(j),:);
        d(1) = d(1) - param.Lx * round(d(1) / param.Lx);
        d(2) = d(2) - param.Ly * round(d(2) / param.Ly);
        diagnostics.minEdgeLength = min(diagnostics.minEdgeLength, norm(d));
    end
end
if isinf(diagnostics.minEdgeLength)
    diagnostics.minEdgeLength = NaN;
end

end
