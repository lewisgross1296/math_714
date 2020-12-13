%% Lewis Gross Math 714 Project
% Implementing Chebyshev Spectral method for solving the 1D Poisson
% Equation

clear ; clc; close all;

%% visualize chebyshev grid and compute DN
N = 50;
[DN,cheb_grid] = cheb(N); 
y =  sqrt(1-cheb_grid.^2);
figure(1);plot(cheb_grid,y,'b-o')
axis([-1.1 1.1 0  1.1])
for j = 1:N
    line([cheb_grid(j),cheb_grid(j)], [0,y(j)]);
end

% forcing functios and analytical solutions, use anonymous style definition
f1 = @(x) 1 - x.^2;
u1 = @(x) (6*x.^2 - x.^4 - 5)/12 ;

f2 = @(x) cos(pi*x) ;
u2 = @(x) -(cos(pi*x) + 1) / pi.^2;

lambda = 4;
f3 = @(x) lambda*x;
u3 = @(x) lambda/6*(x.^3-x);

a = 0.6;
f4 = @(x) heaviside(x+a) - heaviside(x-a) ;
u4 = @(x) (0.5*x.^2 + a.*x + a.^2 ).* heaviside(x+a) - (0.5*x.^2 - a.*x + a.^2 ).* heaviside(x-a) -a*(x+1);

idx=1;
figure(idx);plot(cheb_grid,f1(cheb_grid),'b-o',cheb_grid,f2(cheb_grid),'r-o',...
     cheb_grid,f3(cheb_grid),'k-o',cheb_grid,f4(cheb_grid),'g-o')
xlabel('x');
ylabel('f(x)')
title('forcing functions f(x) on the chebyshev grid')
legend('f1=1-x^2','f2=cos(\pi x)','f3=\lambda x','f4=H(x+a)-H(x-a)','Location','southeast')
idx = idx + 1;

figure(idx);plot(cheb_grid,u1(cheb_grid),'b-o',cheb_grid,u2(cheb_grid),'r-o',...
     cheb_grid,u3(cheb_grid),'k-o',cheb_grid,u4(cheb_grid),'g-o')
xlabel('x');
ylabel('u(x)')
title('analytical solutions on the chebyshev grid')
legend('u1','u2','u3','u4','Location','northwest')
idx = idx + 1;

%% numerical solution, first step slash and compare to analytical

v_slash1 = zeros(N+1,1);
v_slash2 = zeros(N+1,1);
v_slash3 = zeros(N+1,1);
v_slash4 = zeros(N+1,1);

f1_trunc = f1(cheb_grid(2:N));
f2_trunc = f2(cheb_grid(2:N));
f3_trunc = f3(cheb_grid(2:N));
f4_trunc = f4(cheb_grid(2:N));

% make Dn tilde
DN_sq=DN^2;
DN_sq_tilde = DN_sq(2:N,2:N);
v_slash1(2:N) =DN_sq_tilde\f1_trunc;
v_slash2(2:N) =DN_sq_tilde\f2_trunc;
v_slash3(2:N) =DN_sq_tilde\f3_trunc;
v_slash4(2:N) =DN_sq_tilde\f4_trunc;

% numerical solutions
tol = 1e-8;

% WORKING
v_gmres1 = zeros(N+1,1);
v = gmres(@(x) spec_second_deriv_4gmres(x,N),f1_trunc,N-3,tol); % solves cheb mattrix sys
v_gmres1(2:end-1) = v;

figure(idx);plot(cheb_grid,v_gmres1,'bo',cheb_grid,u1(cheb_grid),'r-');
title('Comparing GMRES with Anonymous Function Declaration to Analytical Solution')
legend('GMRES','u1(x)')
idx = idx + 1;

% WORKING
v_gmres2 = zeros(N+1,1);
v = gmres(@(x) spec_second_deriv_4gmres(x,N),f2_trunc,N-3,tol); % solves cheb mattrix sys
v_gmres2(2:end-1) = v;

% WORKING
figure(idx);plot(cheb_grid,v_gmres2,'bo',cheb_grid,u2(cheb_grid),'r-');
title('Comparing GMRES with Anonymous Function Declaration to Analytical Solution')
legend('GMRES','u2(x)')
idx = idx + 1;

% SEEMS TO WORK 
v_gmres3 = zeros(N+1,1);
v = gmres(@(x) spec_second_deriv_4gmres(x,N),f3_trunc,N-3,tol); % solves cheb mattrix sys
v_gmres3(2:end-1) = v;

figure(idx);plot(cheb_grid,v_gmres3,'bo',cheb_grid,u3(cheb_grid),'r-');
title('Comparing GMRES with Anonymous Function Declaration to Analytical Solution')
legend('GMRES','u3(x)')
idx = idx + 1;

v_gmres4 = zeros(N+1,1);
v = gmres(@(x) spec_second_deriv_4gmres(x,N),f4_trunc,N-3,tol); % solves cheb mattrix sys
v_gmres4(2:end-1) = v;
figure(idx);plot(cheb_grid,v_gmres4,'bo',cheb_grid,u4(cheb_grid),'r-');
title('Comparing GMRES with Anonymous Function Declaration to Analytical Solution')
legend('GMRES','u4(x)')
idx = idx + 1;


% show convergence of the analytical solution vs N, determined from
% reguluarity, see trefethen

% how many grid points for a fixed error of GMRES

% show complexity of FFT to compute derivative is NlogN


%% convergence study of GMRES to solve
% show convergence of GMRES, how many iterations required for a fixed error
% fix the number of inner iterations to do this properly
% possible, how costly is GMRES? tic and toc, need to do some large N
% values to show this, i think this is O(N^2)??


%% study how a affects convergence for heaviside 

%% condition number study
NN = [3:1:100];
conds_DN=zeros(size(NN));
conds_DNsq = zeros(size(NN));
conds_DN_til=zeros(size(NN));
conds_DNsq_til = zeros(size(NN));

for n = 1:length(NN)
    [DN,~] = cheb(NN(n)); 
    DN_sq = DN^2;
    conds_DN(n) = cond(DN);
    conds_DNsq(n) = cond(DN_sq);
    conds_DN_til(n) = cond(DN(2:end-1,2:end-1));
    conds_DNsq_til(n) = cond(DN_sq(2:end-1,2:end-1));
end

conds_DN= log(conds_DN);
conds_DNsq = log(conds_DNsq);
conds_DN_til= log(conds_DN_til);
conds_DNsq_til = log(conds_DNsq_til);


figure(idx);plot(NN,conds_DN,'b-o',NN,conds_DNsq,'r-o',NN,conds_DN_til,'k-o',NN,conds_DNsq_til,'m-o')
xlabel('N')
ylabel('condition number')
legend('DN','$DN^{2}$','$\tilde{DN}$','$\tilde{DN}^{2}$','Interpreter','latex')
