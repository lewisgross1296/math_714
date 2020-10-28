%% Lewis Gross Math 714 HW 2 
% This script solves the Cartesian 2D wave equation with homogeneous
% Dirichlet bounday conditions

clear; clc;
% Set up mesh
Ms = [60, 100] ; % number grid points in the x or y direction
hs = 1./(Ms-1);
dts = 0.5*hs;
T_max=10;
for k = 1:length(Ms)
    M = Ms(k) ; h=hs(k); times = [0:dts(k):T_max]; N=length(times) ;
    x = [0:h:1] ;  y = [0:h:1] ;
    [X,Y] = meshgrid(x,y);
    u = zeros(M,M,N) ; % initialize x,y,t mesh
    % propagate initial velocity to first time
    for i = 2:M-1 % loop over x position
        for j = 2:M-1 % loop over y position
            % use IC to find u(i,j,t2), recall t1=0, t2=dt, ...
            % loop to compute the rest of the times
            u(i,j,2) = dts(k)*veloc((i-1)*h)*veloc((j-1)*h) ; % + u(i,j,1) if there was some IC for u
        end
    end
    % compute for the rest of the times
    for n = 3:N
        for i = 2:M-1 % loop over x position
            for j = 2:M-1 % loop over y position
                u(i,j,n) = (dts(k)./h).^2 *( u(i+1,j,n-1) + u(i-1,j,n-1) + u(i,j+1,n-1) + ... 
                       u(i,j-1,n-1)+ (2*(h./dts(k)).^2-4)*u(i,j,n-1) ) - u(i,j,n-2) ;
            end
        end
    end
    
    %plotting
    maxZ = max(max(abs(u(:,:))));
    for n=1:N
        s=surf(X,Y,u(:,:,n));
        zlim([-maxZ,maxZ]);
        pause(0.1);
    end
end
 


