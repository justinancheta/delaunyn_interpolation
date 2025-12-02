clear; close all; clc;
tic
X = 2*gallery('uniformdata',[2000 6],0)-1;
Y = sum(X.^2,2);
d = -0.8:0.4:0.8;
[a0,b0,c0,y0,x0,z0] = ndgrid(d,d,d,d,d,d);
XI = [a0(:) b0(:) c0(:) x0(:) y0(:) z0(:)];
YI = griddatan(X,Y,XI);
t(1,1) = toc;

Tri = delaunayn(X);
Neighbors = build_neighbors(Tri);


tic
for ii = size(XI,1):-1:1
    % Get the test point 
    q = XI(ii,:); 
    
    % Find the nearest simplice
    for jj = size(X,1):-1:1
        r(jj,1) = sqrt( sum( (q - X(jj)).^2, 2) );
    end
    ind = find(r == min(r),1);
    
    simplexID = jump_and_walk_nd(X, Tri, Neighbors, q, ind);

    if simplexID ~= 0
        lambda = barycentric_nd(q, X(Tri(simplexID,:),:) );
        YI2(ii,1) = lambda * Y(Tri(simplexID,:));
    else
        YI2(ii,1) = 0;
    end
    
end
t(2,1) = toc;
plot(YI-YI2)