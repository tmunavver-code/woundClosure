function [longAxisDir, shortAxisDir, centroid] = getCellPrincipalAxis(vertexcoords, vertices)
%GETCELLPRINCIPALAXIS Long/short axis direction of a cell polygon via the
%second-moment (shape) tensor of its vertices about the centroid.
%
% vertexcoords must already be periodicity-corrected for this cell (see
% modifyVerticesForPeriodicity.m) before calling this.
%
% Returns unit vectors longAxisDir (larger eigenvalue -- the cell's
% elongation direction) and shortAxisDir (smaller eigenvalue, orthogonal
% to longAxisDir), plus the centroid used.

pts = vertexcoords(vertices, :);
centroid = mean(pts, 1);
d = pts - centroid;

% Second-moment (shape) tensor
M = (d' * d) / size(d,1);

[V, D] = eig(M);
[~, order] = sort(diag(D), 'descend');
V = V(:, order);

longAxisDir = V(:,1)' / norm(V(:,1));
shortAxisDir = V(:,2)' / norm(V(:,2));

end
