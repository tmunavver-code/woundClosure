function dA = getAreaDerivative(vertexcoords, vertices, currVert, nextVert, prevVert)
    % GETAREADERIVATIVE calculates the derivative of the area of a polygon with respect to the coordinates of currVert.
    % The polygon vertices should be in counter-clockwise order.
    r_next = vertexcoords(nextVert, :);
    r_prev = vertexcoords(prevVert, :);
    
    dA = 0.5 * [r_next(2) - r_prev(2), r_prev(1) - r_next(1)];
end
