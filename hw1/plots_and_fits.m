% this script will generate all the plots without having to wait for the
% solver to re-run for each input (which takes a while)
load('working_error.mat');

figure(1);p1=surf(X,Y,Z);
set(p1,'edgecolor','none')
xlabel('x')
ylabel('y')
zlabel('analytical')

figure(2);p2=surf(X,Y,u');
set(p2,'edgecolor','none')
xlabel('x')
ylabel('y')
zlabel('numerical')

figure(3);p3=surf(X,Y,abs(u'-Z));
set(p3,'edgecolor','none')
xlabel('x')
ylabel('y')
zlabel('abs dif')

figure(4);p4=plot(log(hN),log(errors),'bo');
xlabel('log of hx')
ylabel('log of max error')


coefs_x = polyfit(log(hN),log(errors),1)
coefs_y = polyfit(log(hM),log(errors),1)