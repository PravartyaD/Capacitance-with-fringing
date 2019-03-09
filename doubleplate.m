n = 10;
d = 1;
fun = @(x,y) 1./sqrt(x.^2+y.^2);
self_pot = quad2d(fun, -1./(2*n), 1./(2*n), -1./(2*n), 1./(2*n));
v = zeros(2*n^2);
for k = 0:1:n^2-1
    i_0 = floor(k./n);
    j_0 = mod(k,n);
    x_0 = (j_0./n) + (1./(2*n));
    y_0 = (i_0./n) + (1./(2*n));
    for k1 = 0:1:n^2-1
        i = floor(k1./n);
        j = mod(k1,n);
        if v(k+1,k1+1)==0
            if k == k1
                v(k+1,k1+1) = self_pot;
            else
                func = @(x,y) 1./sqrt((x-x_0).^2+(y-y_0).^2);
                v(k+1,k1+1) = quad2d(func, j./n, (j+1)./n, i./n, (i+1)./n);
                v(k1+1,k+1) = v(k+1,k1+1);
            end
        end
    end
    for k1 = 0:1:n^2-1
        i = floor(k1./n);
        j = mod(k1,n);
        func = @(x,y) 1./sqrt((x-x_0).^2+(y-y_0).^2+d^2);
        v(k+1,k1+1+n^2) = quad2d(func, j./n, (j+1)./n, i./n, (i+1)./n);
    end
end
v = [v(1:n^2,:);v(1:n^2,1+n^2:2*n^2) v(1:n^2,1:n^2)];
v = (9e9/n^2).*v;
V = [0.5*ones(n.^2, 1); -0.5.*ones(n.^2, 1)];
Q = v\V;
Capacitance = sum(Q(1:n^2))
Q1=reshape(Q(1:n^2),[n,n]);
[X,Y]=meshgrid(1/n:1/n:1,1/n:1/n:1);
figure(1)
surf(X,Y,Q1);
hold on;
Q2=reshape(Q(n^2+1:2*n^2),[n,n]);
surf(X,Y,Q2);
hold off;
[X Y] = meshgrid(-0.5:.1:1.5, 0:0.1:d);
v = 0;
for k = 1:1:n^2
    i = floor(k./n);
    j = mod(k,n);
    x = (j./n) + (1./(2*n));
    y = (i./n) + (1./(2*n));
    v = v + Q(k)./sqrt((X-x).^2+(y-0.5).^2+Y.^2) + Q(k+n^2)./sqrt((X-x).^2+(y-0.5).^2+(Y-d).^2);   %Potential Centre Section
end
[px py] = gradient(v);
figure(2)
quiver(X,Y,px,py);
hold off;
