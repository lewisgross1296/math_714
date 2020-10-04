function u = Jacobi_2D_lap_mixed_BCs(ep,M,N)
% this function uses Jacobi iteration to find the solution to 
% the Poisson equation with a dirichlet BC or a one variable boundary 
% distribution for the x boundaries and Neumann conditions for the y
% boundaries

% u is a matrix for the solution at each grid point
% Make initial matrix for Jacobi
u = zeros(M,N);
error=ep+1;
% BCs
% zero at x=1, last column already zero
% first column of u corresponds to x=0 where we have the
% cos(2*pi*y) function, from  the grid y= (j-1)hy
for j = 1:M
    u(j,1) = cos(2*pi*(j-1)*hy) ;
end


while(error > ep)
    %do iteration
end

end
