%% Lewis Gross Math 714 HW 2 
% This script solves the Cartesian 2D wave equation with homogeneous
% Dirichlet bounday conditions

clear; clc;
% Set up mesh
Ms = [60, 120, 180, 240, 360] ; % number grid points in the x or y direction
hs = 1./(Ms-1);

M_fine = 500;
h_fine = 1/(M_fine -1 );
dt = 0.5*h_fine;  % use a dt that satisfies the CFL for the finest grid, will also satisfy for the coarser
T_max=10 ;
times = [0:dt:T_max] ; N=length(times) ;
errors = zeros(length(Ms),N) ;
%fine grid for comparison
x_fine = [0:h_fine:1] ;  y_fine = [0:h_fine:1] ;
[X_fine,Y_fine] = meshgrid(x_fine,y_fine);
u_fine = zeros(M_fine,M_fine,N) ; % initialize x,y,t mesh
% propagate initial velocity to first time
for i = 2:M_fine-1 % loop over x position
    for j = 2:M_fine-1 % loop over y position
        % use IC to find u(i,j,t2), recall t1=0, t2=dt, ...
        % loop to compute the rest of the times
        u_fine(i,j,2) = dt*veloc((i-1)*h_fine)*veloc((j-1)*h_fine) ; % + u(i,j,1) if there was some IC for u
    end
end
    
% compute for the rest of the times
for n = 3:N
    for i = 2:M_fine-1 % loop over x position
        for j = 2:M_fine-1 % loop over y position
            u_fine(i,j,n) = (dt./h_fine).^2 *( u_fine(i+1,j,n-1) + u_fine(i-1,j,n-1) + u_fine(i,j+1,n-1) + ... 
                   u_fine(i,j-1,n-1)+ (2*(h_fine./dt).^2-4)*u_fine(i,j,n-1) ) - u_fine(i,j,n-2) ;
        end
    end
end

% coarse grid solutions
for k = 1:length(Ms)
    M = Ms(k) ; h=hs(k); times = [0:dt:T_max]; N=length(times) ;
    x = [0:h:1] ;  y = [0:h:1] ;
    [X,Y] = meshgrid(x,y);
    u = zeros(M,M,N) ; % initialize x,y,t mesh
    % propagate initial velocity to first time
    for i = 2:M-1 % loop over x position
        for j = 2:M-1 % loop over y position
            % use IC to find u(i,j,t2), recall t1=0, t2=dt, ...
            % loop to compute the rest of the times
            u(i,j,2) = dt*veloc((i-1)*h)*veloc((j-1)*h) ; % + u(i,j,1) if there was some IC for u
        end
    end
    % compute for the rest of the times
    for n = 3:N
        for i = 2:M-1 % loop over x position
            for j = 2:M-1 % loop over y position
                u(i,j,n) = (dt./h).^2 *( u(i+1,j,n-1) + u(i-1,j,n-1) + u(i,j+1,n-1) + ... 
                       u(i,j-1,n-1)+ (2*(h./dt).^2-4)*u(i,j,n-1) ) - u(i,j,n-2) ;
            end
        end
    end
    
    % I HIGHILY RECOMMEND DOING THIS ONCE, IT IS COOL AF
    % uncomment for time evolution plot plotting
%     maxZ = max(max(abs(u(:,:))));
%     for n=1:N
%         s=surf(X,Y,u(:,:,n));
%         zlim([-maxZ,maxZ]);
%         pause(0.05);
%     end

    % error computation, use interp2 to compute the interpolation of the
    % solution on the fine grid to the coarse grid for error analysis
    for n=1:N
        u_fine_to_coarse = interp2(X_fine,Y_fine,u_fine(:,:,n),X,Y);
        errors(k,n) =  max( max(abs(u_fine_to_coarse(:,:) - u(:,:,n)) ) );
    end

end
plot(log(hs),log(errors),'b')
xlabel('log of step size')
ylabel('log of error')