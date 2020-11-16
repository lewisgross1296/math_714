function visualize(u,X,Y)
% generate surf plot for every time (time is the third dimension)
% u is function evaluted at grid X,Y where X,Y are outputs from a meshgrid
% fix picture bounds
    zmax = max(max(abs(u(:,:))));
    sizeu = size(u) ;
    for n=1:sizeu(3) % the number of times in the problem is sizeu(3)
        surf(X,Y,u(:,:,n))
        zlim([-zmax,zmax]);
        pause(0.05)
    end
end