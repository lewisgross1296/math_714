%% Lewis Gross
% This is a script to solve the 2D Poisson equations with the BCS
clear

% Grid and Analytical Solution
a = 0; b = 1 ;
M = 8; N = 8;
hx = (b-a)/(N-1);
hy = (b-a)/(M-1);
x = [a:hx:b] ;
y = [a:hy:b] ;
[X Y] = meshgrid(x,y);
Z = analytical(X,Y) ;
surf(X,Y,Z)


% convergence tolerance for Jacobi
ep = 1e-3; 
u = Jacobi_2D_lap_mixed_BCs(ep,M,N) ;



