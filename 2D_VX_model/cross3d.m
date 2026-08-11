function c = cross3d(a, b)
% CROSS3D Compute cross product of two 3D vectors.
%   Works identically to MATLAB's built-in cross() for 1x3 or 3x1 vectors.

c = [a(2)*b(3) - a(3)*b(2), ...
     a(3)*b(1) - a(1)*b(3), ...
     a(1)*b(2) - a(2)*b(1)];
end
