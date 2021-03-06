function w = chebfft(v)
% Chebyshev differentiation via FFT. Simple, not optimal
% if v is complex, delete the "real" commands
N = length(v) - 1; if N==0, w=0; return, end
x = cos((0:N)'*pi/N);  
ii = 0:N-1;
v=v(:); %ensure v is a column vector
V = [v; flipud(v(2:N))]; % transform x to theta
U = real(fft(V));
W = real(iffft(1i*[ii 0 1-N:-1]'.*U));
w=zeros(N+1,1);
w(2:N) = -W(2:N)./(sqrt(1-x(2:N).^2)); % transform theta to x
w(1) = sum(ii'.^2.*U(ii+1))/N + 0.5*N*U(N+1) ;
w(N+1) = sum((-1).^(ii+1)'.*ii'.^2.*U(ii+1))/N + 0.5*(-1)^(N+1)*N*U(N+1);
end