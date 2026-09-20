function [longAxisDir, shortAxisDir, centroid] = getCellPrincipalAxis(vertexcoords, vertices)
% GETCELLPRINCIPALAXIS Computes principal long and short axis directions and
% centroid of a cell polygon using the second-moment shape tensor.

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
