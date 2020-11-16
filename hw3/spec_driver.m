%% Lewis Gross HW3 Spectral Solver
% inspired from program 20 in Trefethen, Spectral Methods in Matlab

% various N to get convergence as a function of number of grid points
Ns = [400 500 600 650 700 ]; 
% require coarsest grid to be minimum half as many nodes
% put more points closer to the finest grid to sample more in the
% asymptotic region
T = 200;
errors = zeros(size(Ns));

% find grid solution
N_fine = 800;
x_fine = cos(pi*(0:N_fine)/N_fine);
y_fine = x_fine';
[X_fine,Y_fine] = meshgrid(x_fine,y_fine);
u_fine = spec_solver(N_fine,T);

for n=1:length(Ns)
    % compute solution for new N
    [u,Xc,Yc] = spec_solver(Ns(n),T);
    % comment in for visualization
    % visualize(u,Xc,Yc);
    
    % compute the error at each time step, keep the error if it is larger
    % than the error at the previous timestep
    for t=1:T
        % project coarse grid (c) to find grid for error analysis at each
        % time in the solution matrix
        u_proj_coarse2fine = interp2(Xc,Yc,u(:,:,t),X_fine,Y_fine);
        error_t = max(max(abs(u_proj_coarse2fine - u_fine(:,:,t)) ));
        if (error_t>errors(n))
            errors(n) = error_t;
        end
    end
end
coef = polyfit(log(Ns),log(errors),1);
plot(log(Ns),log(errors),'-bo')
xlabel('log of N for Cheb Grid','FontSize',13)
ylabel('log of error','FontSize',13)
title('Spectral Method: Showing Fourth Order Convergence Modified Equation Scheme','FontSize',15)



