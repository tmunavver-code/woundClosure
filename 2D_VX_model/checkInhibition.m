function flag = checkInhibition(indtoadd,coordx,coordy,R_inhibition, Lx, Ly)
% flag = 0: No collision, flag = 1: Collision
flag = 0;
for i = 1:indtoadd-1
	flag = checkCollisionPeriodic(indtoadd,i,coordx,coordy,R_inhibition,Lx,Ly);
	if flag == 1 % Collision detected try again
		break;
	end
end

end

function flag = checkCollisionPeriodic(p1,p2,coordx,coordy,R_inhibition,Lx,Ly)
x1 = coordx(p1); y1 = coordy(p1);
x2 = coordx(p2); y2 = coordy(p2);
dist = getPeriodicDistance(x1,y1,x2,y2, Lx, Ly);
if dist <= 2*R_inhibition
	flag = 1;
else
	flag = 0;
end
end

function dist = getPeriodicDistance(x1,y1,x2,y2,Lx,Ly)
dx = abs(x1 - x2);
if (dx > Lx/2)
	dx = Lx - dx;
end

dy = abs(y1 - y2);
if (dy > Ly/2)
	dy = Ly - dy;
end

dist = sqrt(dx^2 + dy^2);

end