function A = getTriangleArea(vecp1,vecp2,vecp3)

a = vecp2 - vecp1;
b = vecp3 - vecp1;
A = 1/2*norm(cross3d(a,b));