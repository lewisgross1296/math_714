function v_xx = spec_second_deriv(v)
    % accepts function, comptues second derivative
    N = length(v) - 1; 
    if N==0, w=0; return, end
    x = cos((0:N)'*pi/N); 
    v=v(:);
    % transform x to theta
    V = [v; flipud(v(2:N))]; % last element of v is not included
    V_hat = real(fft(V));
    W1 = real(ifft(1i*[0:N-1 0 1-N:-1]'.*V_hat)); % fourier space derivative, ifft to get dtheta
    W2 = real(ifft(-[0:N-1 0 1-N:-1]'.^2.*V_hat)); % fourier spaces second derivative, ifft to get dtheta
    v_xx = zeros(N+1,1);
    v_xx(2:N) = -x(2:N)./(1-x(2:N).^2).^(1.5).*W1(2:N) +  1./(1-x(2:N).^2).*W2(2:N);
    % to do, v_xx(1),v_xx(N+1), do i need these?
    end