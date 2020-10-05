%% Lewis Gross
% This is a script to solve the 2D Poisson equations with the BCS
% with the modified condition for C)(f) 
clear

% list of grid sizes to compute
% x mesh point number
NN = [6,10,22,40,52];
hN = 1./(NN-ones(1,length(NN)));
% y mesh point number
MM = [4,8,20,38,50];
hM = 1./(MM-ones(1,length(MM)));
a = 0; b = 1 ;
errors=zeros(length(NN),1);
for k = 1:length(NN)
    % Grid and Analytical Solution
    N = NN(k); M =MM(k);
    hx = hN(k);
    hy = hM(k);
    x = [ a : hx : b ] ;  % x1, x2, ..., xi, ... , xN-1, xN
    y = [ a : hy : b ] ;  % y1, x2, ..., xi, ... , xN-1, xN
    [X , Y] = meshgrid(x,y);
    Z = analytical(X,Y) ; 
  
    % convergence tolerance for Jacobi, much lower since we're expecting a
    % bad result, no need to have the iterations go for as long as before
    ep = 1e-5; 
    % call solver
    u = Jacobi_bad_BC(ep,N,M,hx,hy);
    errors(k) = max( max(abs(u'-Z)) );  %compute max error among all grid points
end
figure(1);plot(log(NN),log(errors),'ro')
xlabel('log of N grid points')
ylabel('log of error')

figure(2);p2=surf(X,Y,u')
set(p2,'edgecolor','none')
xlabel('x')
ylabel('y')
zlabel('numerical')

figure(3);p3=surf(X,Y,abs(u'-Z))
set(p3,'edgecolor','none')
xlabel('x')
ylabel('y')
zlabel('abs dif')