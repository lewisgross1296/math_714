%% Lewis Gross Math 714 Project
% Implementing Chebyshev Spectral method for solving the 1D Poisson
% Equation

clear ; clc;

% visualize chebyshev grid and compute DN
N = 16;
[D,cheb_grid] = cheb(N); 
y =  sqrt(1-cheb_grid.^2);
figure(1);plot(cheb_grid,y,'b-o')
axis([-1.1 1.1 0  1.1])
for j = 1:N
    line([cheb_grid(j),cheb_grid(j)], [0,y(j)]);
end


% forcing functinos and analytical solutions, use anonymous style definition

f1 = @(x) 1 - x.^2;
u1 = @(x) (6*x.^2 - x.^4 - 5)/12 ;

f2 = @(x) cos(pi*x/2) ;
u2 = @(x) -(cos(pi*x) + 1) / pi.^2;

lambda = 4;
f3 = @(x) lambda*x;
u3 = @(x) lambda/6*(x.^3-x);

a = 0.6;
f4 = @(x) heaviside(x+a) - heaviside(x-a) ;
u4 = @(x) (0.5*x.^2 + a.*x + a.^2 ).* heaviside(x+a) - (0.5*x.^2 - a.*x + a.^2 ).* heaviside(x-a);

figure(2);plot(cheb_grid,f1(cheb_grid),'b-o',cheb_grid,f2(cheb_grid),'r-o',...
     cheb_grid,f3(cheb_grid),'k-o',cheb_grid,f4(cheb_grid),'g-o')
xlabel('x');
ylabel('f(x)')
title('forcing functions f(x) on the chebyshev grid')
legend('f1','f2','f3','f4')

figure(3);plot(cheb_grid,u1(cheb_grid),'b-o',cheb_grid,u2(cheb_grid),'r-o',...
     cheb_grid,u3(cheb_grid),'k-o',cheb_grid,u4(cheb_grid),'g-o')
xlabel('x');
ylabel('u(x)')
title('analytical solutions on the chebyshev grid')
legend('u1','u2','u3','u4')

% numerical solution, first step, convergence vs N

% make Dn tilde
% vectorize f
% slash it up


% numerical solution, 
