function [D,cheb_grid] = cheb(N)
% From Trefethen, SIAM 200 chem.m (chapter 6 page 54)
if (N==0)
    D=0; cheb_grid = 1;
    return
else
    cheb_grid = cos([0:1:N]*pi/N)'; % represented as a column vector for the following
    c = [2 ; ones(N-1,1) ; 2] .*(-1).^[0:N]'; 
    X = repmat(cheb_grid,1,N+1); % copies cheb_grid column N times to make a (N+1)x(N+1) matrix
    dX = X - X';
    D = (c*(1./c)')./(dX + eye(N+1));
    D = D - diag(sum(D')); % higher accuracy way to compute the diagonal
end

end
