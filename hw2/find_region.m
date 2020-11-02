x = [-50:1:50];
y = [-50:1:50];
[real , imag] = meshgrid(x,y);
alpha = real+sqrt(-1)*imag;


rho_plus = 1 + alpha/2 + sqrt((1 + alpha/2).^2-1);
rho_minus = 1 - alpha/2 - sqrt((1 + alpha/2).^2-1);

mag_plus = abs(rho_plus);
mag_minus = abs(rho_minus);
onesurf = ones(size(alpha));
surf(real,imag,onesurf);




