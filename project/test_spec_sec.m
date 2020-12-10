% testing spec second deriv
N = 100;
x = cos((0:N)'*pi/N);
v1 = sin(x);
v_xx1 = spec_second_deriv(v1);
figure(1);plot(x,v1,'bo-',x,v_xx1,'ro-')

v2 = x.^2;
v_xx2 = spec_second_deriv(v2);
figure(2);plot(x,v2,'bo-',x,v_xx2,'ro-')

v3 = cos(x);
v_xx3 = spec_second_deriv(v3);
figure(3);plot(x,v3,'bo-',x,v_xx3,'ro-')

v4= exp(x);
v_xx4 = spec_second_deriv(v4);
figure(4);plot(x,v4,'bo-',x,v_xx4,'ro-')

v5= 4*x.*sin(4*x);
v_xx5 = spec_second_deriv(v5);
figure(5);plot(x,v5,'bo-',x,v_xx5,'ro-')