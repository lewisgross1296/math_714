function [u , X, Y] = spec_solver(N,T)
% This function accepts a number of grid points for the Chebyshev grid N
% and a number of time steps T
% It solves the the wave equation with the initial velocity profile
% specified in veloc.m

% It uses a correction to the wave equation to make a fourth order in time
% method for u(dt) for the second order in time method approximation to the
% second time derivative

x = cos(pi*(0:N)/N);
y = x';
dt = 6/N^2;
[X,Y] = meshgrid(x,y);

plotgap = round((1/3)/dt);
dt = (1/3)/plotgap;
tfin = T*dt ; % time simulation ends

u = zeros(N+1,N+1,T); % 3d matrix
ut0 =veloc_domain_shift(x).*veloc_domain_shift(y); 
% initial velocity requires shifted domain, as our problem is [0,1]^2 but
% the spectral method requires [-1,1]^2
% work was done to transform the velocity funciton to behave correctly for
% the shifted grid (as compared to the FD driver used in HW2

% fourth order initialization at t=dt
u(:,:,2) = u(:,:,1) + dt*ut0 + dt^2/2*laplacian(u(:,:,1),x,y,N) ...  
       + dt^3/6*laplacian(ut0,x,y,N) ;
   
% for the rest of the times, compute the solution
for n = 3:T
    lap = laplacian(u(:,:,n-1),x,y,N);
    bilap = laplacian(lap,x,y,N);
    u(:,:,n) = 2*u(:,:,n-1) - u(:,:,n-2) + dt^2*(lap) + dt^4/12*bilap;
end


end
