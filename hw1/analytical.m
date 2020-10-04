function [u] = analytical(x,y)
% analytical solution to laplaces equation with BCS
% u(0,y) = cos(2piy)
% u(1,y) = 0
% u_y(x,0)  = 0
% u_y(x,1)  = 0

u = cos(2*pi*y).*(cosh(2*pi.*x)- coth(2*pi).*sinh(2*pi.*x));
end