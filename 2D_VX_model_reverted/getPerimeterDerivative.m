function dPeri = getPerimeterDerivative(vertexcoords, currVert, nextVert, prevVert)
    % GETPERIMETERDERIVATIVE calculates the derivative of the perimeter with respect to the coordinates of currVert.
    r_curr = vertexcoords(currVert, :);
    r_next = vertexcoords(nextVert, :);
    r_prev = vertexcoords(prevVert, :);
    
    vec1 = r_curr - r_prev;
    vec2 = r_curr - r_next;
    
    L1 = norm(vec1);
    L2 = norm(vec2);
    
    dPeri = zeros(1, 2);
    if L1 > 0
        dPeri = dPeri + vec1 / L1;
    end
    if L2 > 0
        dPeri = dPeri + vec2 / L2;
    end
end
