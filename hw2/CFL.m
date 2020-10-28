%% script to plot CFL condition
% still some things to work out

clear; clc;
kh = linspace(0,45,1000);
k = [1 , 20, 50, 100 , 200, 1000];
z = zeros(size(kh));
dt = 0.01;
for it=1:length(k)
    h = kh/k(it);
    upper = 2 * h.^2./dt.^2;
    CFL_plus = @(kh) 2 - cos(kh) + sqrt( (cos(kh) + h.^2/dt.^2 -2 ).^2 - 1 ) ;
    CFL_minus = @(kh) 2 - cos(kh) - sqrt( (cos(kh) + h.^2/dt.^2 -2 ).^2 - 1 ) ;
    figure(it);plot(kh,z,'k-',kh,CFL_plus(kh),'r',kh,CFL_minus(kh),'b',kh, upper, 'k-.')
    title(['CFL plot of bounds as a function of kh for \Deltat=',num2str(dt),...
        ' and a value of k=',num2str(k(it))])
    xlabel('kh')
    ylabel('CFL function and bounds')
    legend('lower bound','plus root','minus root','upper bound')
end

dt = 0.1;
for it=1:length(k)
    h = kh/k(it);
    upper = 2 * h.^2./dt.^2;
    CFL_plus = @(kh) 2 - cos(kh) + sqrt( (cos(kh) + h.^2/dt.^2 -2 ).^2 - 1 ) ;
    CFL_minus = @(kh) 2 - cos(kh) - sqrt( (cos(kh) + h.^2/dt.^2 -2 ).^2 - 1 ) ;
    figure(it);plot(kh,z,'k-',kh,CFL_plus(kh),'r',kh,CFL_minus(kh),'b',kh, upper, 'k-.')
    title(['CFL plot of bounds as a function of kh for \Deltat=',num2str(dt),...
        ' and a value of k=',num2str(k(it))])
    xlabel('kh')
    ylabel('CFL function and bounds')
    legend('lower bound','plus root','minus root','upper bound')
end

dt = 0.5;
for it=1:length(k)
    h = kh/k(it);
    upper = 2 * h.^2./dt.^2;
    CFL_plus = @(kh) 2 - cos(kh) + sqrt( (cos(kh) + h.^2/dt.^2 -2 ).^2 - 1 ) ;
    CFL_minus = @(kh) 2 - cos(kh) - sqrt( (cos(kh) + h.^2/dt.^2 -2 ).^2 - 1 ) ;
    figure(it);plot(kh,z,'k-',kh,CFL_plus(kh),'r',kh,CFL_minus(kh),'b',kh, upper, 'k-.')
    title(['CFL plot of bounds as a function of kh for \Deltat=',num2str(dt),...
        ' and a value of k=',num2str(k(it))])
    xlabel('kh')
    ylabel('CFL function and bounds')
    legend('lower bound','plus root','minus root','upper bound')
end