function woundArea = getWoundArea(celldata, param)
%GETWOUNDAREA Total wound area, by conservation of the periodic box:
%
%   sum(cell areas) + wound area = Lx * Ly
%
% The tissue tiles the whole periodic domain apart from the wound hole,
% so the wound area is just whatever box area the cells do not cover.
%
% This replaces measuring the wound by tracing its boundary loop
% (getOrderedWoundLoop.m), which assumes the margin is a single SIMPLE
% closed loop -- every vertex of degree exactly 2. That assumption breaks
% legitimately as a wound closes: the boundary can narrow to a waist and
% pinch at a vertex, leaving that vertex with degree 4 and the wound in
% two lobes. Loop tracing then reports "invalid" and the area trace goes
% NaN for the rest of the run, even though the mesh is perfectly
% consistent (verified: verttocell agreed with connec, and the pinch
% vertex had even degree, at the step where this first occurred).
%
% Measuring by area conservation has no topology assumptions at all, so
% it keeps working through pinches, multiple lobes, and final closure.

woundArea = param.Lx * param.Ly - sum(celldata.A);

if woundArea < 0
    woundArea = 0; % numerically closed
end

end
