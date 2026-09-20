function [coordinates, connectivity, verttocell, ok] = performWoundSeal(v1, v2, coordinates, connectivity, verttocell)
% PERFORMWOUNDSEAL Merges two opposing wound-margin vertices (v1, v2) into
% their midpoint, retiring boundary edges into shared internal edges.
% Returns ok=false if the merge produces a degenerate cell.

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
