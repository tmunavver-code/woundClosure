function [tilex, tiley] = getTile(coordx,coordy,Lx,Ly)
% Generate 9X9 tile
N = length(coordx);
tilex = zeros(9*N,1);
tiley = zeros(9*N,1);

% Tile 1
tile = 1;
tilex(N*(tile-1)+1:N*tile) = coordx;
tiley(N*(tile-1)+1:N*tile) = coordy;

% Tile 2
tile = 2;
tilex(N*(tile-1)+1:N*tile) = coordx - Lx;
tiley(N*(tile-1)+1:N*tile) = coordy + Ly;

% Tile 3
tile = 3;
tilex(N*(tile-1)+1:N*tile) = coordx + Lx;
tiley(N*(tile-1)+1:N*tile) = coordy + Ly;

% Tile 4
tile = 4;
tilex(N*(tile-1)+1:N*tile) = coordx + Lx;
tiley(N*(tile-1)+1:N*tile) = coordy - Ly;

% Tile 5
tile = 5;
tilex(N*(tile-1)+1:N*tile) = coordx - Lx;
tiley(N*(tile-1)+1:N*tile) = coordy - Ly;

% Tile 6
tile = 6;
tilex(N*(tile-1)+1:N*tile) = coordx;
tiley(N*(tile-1)+1:N*tile) = coordy + Ly;

% Tile 7
tile = 7;
tilex(N*(tile-1)+1:N*tile) = coordx + Lx;
tiley(N*(tile-1)+1:N*tile) = coordy;

% Tile 8
tile = 8;
tilex(N*(tile-1)+1:N*tile) = coordx;
tiley(N*(tile-1)+1:N*tile) = coordy - Ly;

% Tile 9
tile = 9;
tilex(N*(tile-1)+1:N*tile) = coordx - Lx;
tiley(N*(tile-1)+1:N*tile) = coordy;

end



