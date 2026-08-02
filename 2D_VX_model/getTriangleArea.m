function A = getTriangleArea(p1, p2, p3)
    % GETTRIANGLEAREA computes the area of a triangle defined by p1, p2, p3.
    % p1, p2, p3 are expected to be 1x3 vectors (e.g., [x, y, 0]).
    v1 = p2 - p1;
    v2 = p3 - p1;
    c = cross(v1, v2);
    A = 0.5 * norm(c);
end
