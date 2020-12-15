%% Lewis Gross Math 714 Project
% Implementing Chebyshev Spectral method for solving the 1D Poisson
% Equation

clear ; clc; close all;

%% visualize chebyshev grid and compute DN
N = 50;
[DN,cheb_grid] = cheb(N); 
y =  sqrt(1-cheb_grid.^2);
idx = 1;
figure(idx);plot(cheb_grid,y,'b-o')
title(['Chebyshev Grid for N=',num2str(N)])
xlabel('Grid Locations')
ylabel('Projection onto the Unit Cirlce')
idx = idx  + 1;

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

%% slash solution compared to analytical
close all;
% % switch for each f1, f2, or f3
% % Ns = [128 145 152 161 184 190]; % good for f1
% % Ns = [128 142 149 169 178 197]; % good for f2
% Ns = [120 144 147 165 168 178 189 197]; % good for f3
% errors1 = zeros(size(Ns)); 
% errors2 = zeros(size(Ns));
% errors3 = zeros(size(Ns));
% errors4 = zeros(size(Ns));
% 
% for n=1:length(Ns)
%     
%     N = Ns(n);
%     [DN,cheb_grid] = cheb(N); 
%     v_slash1 = zeros(N+1,1);
%     v_slash2 = zeros(N+1,1);
%     v_slash3 = zeros(N+1,1);
%     v_slash4 = zeros(N+1,1);
% 
%     f1_trunc = f1(cheb_grid(2:N));
%     f2_trunc = f2(cheb_grid(2:N));
%     f3_trunc = f3(cheb_grid(2:N));
%     f4_trunc = f4(cheb_grid(2:N));
% 
%     % make Dn tilde
%     DN_sq=DN^2;
%     DN_sq_tilde = DN_sq(2:N,2:N);
%     v_slash1(2:N) =DN_sq_tilde\f1_trunc;
%     v_slash2(2:N) =DN_sq_tilde\f2_trunc;
%     v_slash3(2:N) =DN_sq_tilde\f3_trunc;
%     v_slash4(2:N) =DN_sq_tilde\f4_trunc;
%     
%     figure(idx);plot(cheb_grid,v_slash1,'bo',cheb_grid,u1(cheb_grid),'r-')
%     title(['Comparing Analytical and Slash Solutions for N=',num2str(N)])
%     legend('slash solution 1','u1(x)')
%     idx = idx + 1;
% 
%     figure(idx);plot(cheb_grid,v_slash2,'bo',cheb_grid,u2(cheb_grid),'r-')
%     title(['Comparing Analytical and Slash Solutions for N=',num2str(N)])
%     legend('slash solution 2','u2(x)')
%     idx = idx + 1;
% 
%     figure(idx);plot(cheb_grid,v_slash3,'bo',cheb_grid,u3(cheb_grid),'r-')
%     title(['Comparing Analytical and Slash Solutions for N=',num2str(N)])
%     legend('slash solution 3','u3(x)')
%     idx = idx + 1;
%     
%     figure(idx);plot(cheb_grid,v_slash4,'bo',cheb_grid,u4(cheb_grid),'r-')
%     title(['Comparing Analytical and Slash Solutions for N=',num2str(N)])
%     legend('slash solution 4','u4(x)')
%     idx = idx + 1;
% 
%     errors1(n) = norm(abs(v_slash1 - u1(cheb_grid)));
%     errors2(n) = norm(abs(v_slash2 - u2(cheb_grid)));
%     errors3(n) = norm(abs(v_slash3 - u3(cheb_grid)));
%     % show via picture that 4 doesn't work, so don't compute errors
%     
% end
% 
% 
% % comment in correct plot to get linear shape 
% % figure(idx); plot(Ns,log10(errors1),'b-o')
% % xlabel('N')
% % ylabel('log(error)')
% % title('O($K^{-N}$) convergence','Interpreter','latex')
% % legend('f1')
% % idx = idx + 1;
% % 
% % figure(idx); plot(Ns,log10(errors2),'r-o')
% % xlabel('N')
% % ylabel('log(error)')
% % title('O($K^{-N}$) convergence','Interpreter','latex')
% % legend('f2')
% % idx = idx + 1;
% 
% figure(idx); plot(Ns,log10(errors3),'m-o') 
% xlabel('N')
% ylabel('log(error)')
% title('O($K^{-N}$) convergence','Interpreter','latex')
% legend('f3')
% idx = idx + 1;

%% Numerical solution without slash, solving system approximately (for case using 3 or 4digis)

Ns = [ 14 15 16 17 18 19 20 22 24];
tol = 1e-3;

% iterations for a fixed error from GMRES
iters_count1 = zeros(size(Ns));
iters_count2 = zeros(size(Ns));
iters_count3 = zeros(size(Ns));

for n = 1:length(Ns)
    N = Ns(n);
    restart = 12; % 20 to avoid useless loss of orthogonality comps
    maxit = 5; % 100 to get to the sol, revisit with N
    [~,cheb_grid] = cheb(N); 
    f1_trunc = f1(cheb_grid(2:N));
    f2_trunc = f2(cheb_grid(2:N));
    f3_trunc = f3(cheb_grid(2:N));

    % f1
    v_gmres1 = zeros(N+1,1);
    [v1 , flag1, relres1, iter1] = gmres(@(x) spec_second_deriv_4gmres(x,N),f1_trunc,restart,tol,maxit); 
    v_gmres1(2:end-1) = v1;
    if(iter1(1)==1)
        iters_count1(n) = iter1(2);
    else
        iters_count1(n) = (iter1(1)-1)*maxit + iter1(2) ;
    end

    % f2
    v_gmres2 = zeros(N+1,1);
    [v2 , flag2, relres2, iter2] = gmres(@(x) spec_second_deriv_4gmres(x,N),f2_trunc,restart,tol,maxit); 
    v_gmres2(2:end-1) = v2;
    idx = idx + 1;
    if(iter2(1)==1)
        iters_count2(n) = iter2(2);
    else
        iters_count2(n) = (iter2(1)-1)*maxit + iter1(2) ;
    end
    % f3
    v_gmres3 = zeros(N+1,1);
    [v3 , flag3, relres3, iter3] = gmres(@(x) spec_second_deriv_4gmres(x,N),f3_trunc,restart,tol,maxit); 
    v_gmres3(2:end-1) = v3;
    if(iter3(1)==1)
        iters_count3(n) = iter3(2);
    else
        iters_count3(n) = (iter3(1)-1)*maxit + iter3(2) ;
    end
    % flag1, flag2, flag3 
    % uncomment to see all flags equaling zero, indicating convergence
end

figure(idx); plot(Ns,iters_count1,'b-o',Ns,iters_count2,'r-+',Ns,iters_count3,'m-s')
xlabel('grid size')
ylabel('iteration count')
legend('f1','f2','f3')
title('study of iterations necessary to reach an error of 1e-3')
idx = idx + 1;

%% condition number study

NN = [20:10:600];
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

figure(idx);semilogy(NN,conds_DN,'b-o',NN,conds_DNsq,'r-o',NN,conds_DN_til,'k-o',NN,conds_DNsq_til,'m-o')
xlabel('N')
ylabel('condition number')
legend('DN','$DN^{2}$','$\tilde{DN}$','$\tilde{DN}^{2}$','Interpreter','latex')
title('Condition Number for Various Forms of Chebyshev Differentiation Matrix')