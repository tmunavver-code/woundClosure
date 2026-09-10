function [coordinates, connectivity, verttocell, ok] = performWoundSeal(v1, v2, coordinates, connectivity, verttocell)
%PERFORMWOUNDSEAL Merge two wound-margin vertices that lie on OPPOSITE
%sides of a collapsing wound, sealing the wound at that point.
%
% This is the step that actually retires cells from the wound margin.
% Without it the wound area can reach zero while the margin still reports
% its full complement of cells and free edges: the hole collapses to a
% degenerate sliver, but nothing ever converts those free edges into
% ordinary shared edges, so no cell ever leaves the boundary. (Measured:
% wound area 0.00% at step 700 with 15 margin cells / 22 free edges still
% reported.)
%
% v1 survives and is moved to the midpoint; v2 is orphaned and every cell
% referencing it is repointed to v1. Free edges on either side that had
% v2 as an endpoint now meet cells from the other side at v1, so they
% become internal edges and drop out of the free-edge set.
%
% Returns ok=false with nothing modified if the merge would produce a
% degenerate (<3 vertex) cell.

ok = false;

cellsAtV2 = verttocell{v2};

% Validate first -- no mutation until every affected cell is known good.
for k = 1:numel(cellsAtV2)
    c = cellsAtV2(k);
    list = connectivity{c};
    list(list == v2) = v1;
    if numel(unique(list)) < 3
        return;
    end
end

% Commit: v1 moves to the midpoint of the two, v2 is absorbed into it.
coordinates(v1, :) = 0.5 * (coordinates(v1, :) + coordinates(v2, :));

for k = 1:numel(cellsAtV2)
    c = cellsAtV2(k);
    list = connectivity{c};
    list(list == v2) = v1;
    connectivity{c} = list;
end

verttocell{v1} = unique([verttocell{v1}, cellsAtV2]);
verttocell{v2} = [];

ok = true;

end
