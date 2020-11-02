%% Lewis Gross Math 714 HW2 Question B
% This script determines an N such that a function and its linear
% interpolant differ by less than 10^(-2)
clear;
tol = 10^(-2);
error = tol  + 1;
N = 10;
% even loop
while(error>tol)
    N=N+1;
    h = 1/(N-1) ;
    xgrid = [0:h:1] ;
    errors = zeros(size(xgrid));
    for j = 1:length(xgrid)-1
        xguess = (xgrid(j)+xgrid(j+1))/2;
        xstar = fzero( @(x) gprime( x , xgrid(j+1) , xgrid(j) , N) , xguess );
        interp_star = veloc(xgrid(j)) + (N-1)*(veloc(xgrid(j+1)) -veloc(xgrid(j))) *( xstar - xgrid(j) ) ;
        errors(j) = veloc(xstar) - interp_star;
    end
    error = max(abs(errors));
end
disp(['N is ',num2str(N)]);
% tol = 10^(-2);
% error = tol  + 1;
% N_odd = 11;
% % even loop
% while(error>tol)
%     h = 1/(N_odd-1) ;
%     xgrid = [0:h:1] ;
%     errors = zeros(size(xgrid));
%     for j = 1:length(xgrid)-1
%         xguess = (xgrid(j)+xgrid(j+1))/2;
%         xstar = fzero( @(x) gprime( x , xgrid(j+1) , xgrid(j) , N_odd) , xguess );
%         interp_star = veloc(xgrid(j)) + (N_odd-1)*(veloc(xgrid(j+1)) -veloc(xgrid(j))) *( xstar - xgrid(j) ) ;
%         errors(j) = veloc(xstar) - interp_star;
%     end
%     error = max(abs(errors));
%     N_odd=N_odd+2;
% end
% 
% disp(['Even N is ',num2str(N),'. Odd N is',num2str(N_odd)])
