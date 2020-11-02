%% Code based off code from Varun Gudibanda
% Find region of Stability for C(b)on HW2
f = @(x,y) 1+(x+i*y)/2 + sqrt((1+(x+i*y)/2).^2 - 1);
g = @(x,y) 1+(x+i*y)/2 - sqrt((1+(x+i*y)/2).^2 - 1);
[X,Y] = meshgrid(-5:0.1:1,-1:0.1:1);
h = ones(size(X));

figure(1)
surf(X, Y, abs(f(X,Y)), 'FaceColor','g', 'FaceAlpha',0.75, 'EdgeColor','none')
hold on
surf(X, Y, abs(g(X,Y)), 'FaceColor','r', 'FaceAlpha',0.75, 'EdgeColor','none')
surf(X, Y, h, 'FaceColor', 'b', 'FaceAlpha', 0.5, 'EdgeColor', 'none')
xlabel('Real')
ylabel('Imag')
hold off