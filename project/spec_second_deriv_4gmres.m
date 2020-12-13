function v_xx = spec_second_deriv_4gmres(v,N)
    % accepts vector that is size of the truncated system
    % comptues second derivative except at the end points of the chebyshev grid
    % N is the size of the cheb grid
    x = cos((0:N)'*pi/N); 
    v=v(:); 
    % transform x to theta
    V = [ 0; v; 0; flipud(v(1:N-1)) ]; % v given by function is missing endpoints v0=vN=0
    V_hat = real(fft(V));
    W1 = real(ifft(1i*[0:N-1 0 1-N:-1]'.*V_hat)); % fourier space derivative, ifft to get dtheta
    W2 = real(ifft(-[0:N-1 0 1-N:-1]'.^2.*V_hat)); % fourier spaces second derivative, ifft to get dtheta
    v_xx = zeros(N-1,1);
    v_xx = -x(2:N)./(1-x(2:N).^2).^(1.5).*W1(2:N) +  1./(1-x(2:N).^2).*W2(2:N);
    % to do, v_xx(1),v_xx(N+1), do i need these?
    end