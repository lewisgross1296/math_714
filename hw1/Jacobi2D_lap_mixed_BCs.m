function u = Jacobi2D_lap_mixed_BCs(ep,N,M,hx,hy)
% this function uses Jacobi iteration to find the solution to 
% the Poisson equation with a dirichlet BC or a one variable boundary 
% distribution for the x boundaries and Neumann conditions for the y
% boundaries

% u is a matrix for the solution at each grid point
% Make initial matrix for Jacobi
u = zeros(N,M);
% rows correspond to x positions
% columns correspond to y positions

% constant for iteration
c = (hx^2*hy^2)/( 2*(hx^2+hy^2) );

% BCs
% zero at x=1, last row already zero
for j = 1:M
    u(1,j) = cos(2*pi*(j-1)*hy) ;
end

u_init = u;
error = ep+1;
while(error > ep)
    unew = u_init;
    % do iteration
        for i = 2:N-1
            % boundary j=1 (y=0), use ghost node
            unew(i,1)= c*( ( u(i-1,1) + u(i+1,1) )/(hx^2) +  (  2*u(i,2))/(hy^2) ) ;
            for j=2:M-1
                unew(i,j) = c*( ( u(i-1,j) + u(i+1,j) )/(hx^2) +  ( u(i,j-1) + u(i,j+1) )/(hy^2) );
            end
            % boundary j=M (y=1), use ghost node
            unew(i,M)= c*( ( u(i-1,M) + u(i+1,M) )/(hx^2) +  (  2*u(i,M-1))/(hy^2) ) ;
        end
        % compute the 2 norm to find the error of the new iteration
        error = norm(unew-u,2);
        % store newest computed solution as the solution to be retturned,
        % if it is converged, this u will be returned 
        u = unew;
end
end