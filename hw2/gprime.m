function g = gprime(x,xhigh,xlow,N)
% x is the position, xhigh is the upper bound for the interval ,xlow is
% the lower bound for the interval, N is the number of grid points in the mesh
    g  = (400-800*x).*exp(-400.*(x-0.5).^2 ) - (N-1).*(veloc(xhigh) - veloc(xlow));
end