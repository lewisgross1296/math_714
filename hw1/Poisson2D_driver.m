%% Lewis Gross
% This is a script to solve the 2D Poisson equations with the BCS
% specified in part C)(a) 
clear

% list of grid sizes to compute
% x mesh point number
NN = [6,22,52,102,202];
hN = 1./(NN-ones(1,length(NN)));
% y mesh point number
MM = [4,20,50,100,200];
hM = 1./(MM-ones(1,length(MM)));
a = 0; b = 1 ;
errors=zeros(length(NN),1);
counts=zeros(length(NN),1);
for k = 1:length(NN)
    % Grid and Analytical Solution
    N = NN(k); M =MM(k);
    hx = hN(k);
    hy = hM(k);
    x = [ a : hx : b ] ;  % x1, x2, ..., xi, ... , xN-1, xN
    y = [ a : hy : b ] ;  % y1, x2, ..., xi, ... , xN-1, xN
    [X , Y] = meshgrid(x,y);
    Z = analytical(X,Y) ; 
  
    % convergence tolerance for Jacobi
    ep = 1e-14; 
    % call solver
    [u, c] = Jacobi2D_lap_mixed_BCs(ep,N,M,hx,hy);
    counts(k) = c;
    errors(k) = max( max(abs(u'-Z)) );  %compute max error among all grid points
end