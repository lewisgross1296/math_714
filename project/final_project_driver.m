%% Lewis Gross Math 714 Project
% Implementing Chebyshev Spectral method for solving the 1D Poisson
% Equation

clear ; clc;

%% visualize chebyshev grid and compute DN
N = 30;
[DN,cheb_grid] = cheb(N); 
y =  sqrt(1-cheb_grid.^2);
figure(1);plot(cheb_grid,y,'b-o')
axis([-1.1 1.1 0  1.1])
for j = 1:N
    line([cheb_grid(j),cheb_grid(j)], [0,y(j)]);
end


%% forcing functinos and analytical solutions, use anonymous style definition

f1 = @(x) 1 - x.^2;
u1 = @(x) (6*x.^2 - x.^4 - 5)/12 ;

f2 = @(x) cos(pi*x/2) ;
u2 = @(x) -(cos(pi*x) + 1) / pi.^2;

lambda = 4;
f3 = @(x) lambda*x;
u3 = @(x) lambda/6*(x.^3-x);

a = 0.6;
f4 = @(x) heaviside(x+a) - heaviside(x-a) ;
u4 = @(x) (0.5*x.^2 + a.*x + a.^2 ).* heaviside(x+a) - (0.5*x.^2 - a.*x + a.^2 ).* heaviside(x-a) -a*(x+1);

figure(2);plot(cheb_grid,f1(cheb_grid),'b-o',cheb_grid,f2(cheb_grid),'r-o',...
     cheb_grid,f3(cheb_grid),'k-o',cheb_grid,f4(cheb_grid),'g-o')
xlabel('x');
ylabel('f(x)')
title('forcing functions f(x) on the chebyshev grid')
legend('f1=1-x^2','f2=cos(\pi x)','f3=\lambda x','f4=H(x+a)-H(x-a)','Location','southeast')

figure(3);plot(cheb_grid,u1(cheb_grid),'b-o',cheb_grid,u2(cheb_grid),'r-o',...
     cheb_grid,u3(cheb_grid),'k-o',cheb_grid,u4(cheb_grid),'g-o')
xlabel('x');
ylabel('u(x)')
title('analytical solutions on the chebyshev grid')
legend('u1','u2','u3','u4','Location','northwest')

%% numerical solution, first step slash and compare to analytical
v = zeros(N+1,1);
f1_trunc = f1(cheb_grid(2:end-1));

figure(4);plot(cheb_grid(2:end-1),f1_trunc,cheb_grid,f1(cheb_grid))
% f2_trunc
% f3_trunc
% f4_trunc
% make Dn tilde squared
DN_sq=DN^2;
DN_sq_tilde = DN_sq(2:end-1,2:end-1);
v(2:end-1) =DN_sq_tilde\f1_trunc;

v(2:end-1) = gmres(DN_sq_tilde,f1_trunc);

% do gmres with function handle that is responsible for computing u''
% see laplacian.m, one go of the loop tho

% show complexity of FFT to compute derivative is NlogN

% show convergence of the analytical solution vs N, determined from
% reguluarity, see trefethen

% show convergence of GMRES, how many iterations required for a fixed error
% fix the number of inner iterations to do this properly
% possible, how costly is GMRES? tic and toc, need to do some large N
% values to show this, i think this is O(N^2)??

% how many grid points for a fixed error of GMRES

figure(5);plot(cheb_grid,v,'bo',cheb_grid,u1(cheb_grid),'r')


%% convergence study of GMRES to solve


%% study how a affects convergence for heaviside 
