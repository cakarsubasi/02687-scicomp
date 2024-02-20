u_exact = @(x, y) sin(4*pi.*(x+y)) + cos(4.*pi.*x.*y);

f = @(x, y) -16*pi^2 * (2*sin(4*pi*(x+y)) + (x.^2 +y.^2)*cos(4*pi.*x.*y));

%% BCs
a = 0;
b = 1;
c = 0;
d = 1;

m = 3;
h = (b-a)/(m+1);

X = linspace(a, b, m+2); % all points
Y = linspace(c, d, m+2);
[X,Y] = meshgrid(X,Y);
X = X';
Y = Y';
Xindices = 2:m+1; % interior points
Yindices = 2:m+1;
Xint = X(Xindices,Yindices);       % interior points
Yint = Y(Xindices,Yindices);

A = poisson5(m);

% solution:
u_solution = u_exact(X, Y);
%%
rhs = f(Xint, Yint);

% adjust the rhs to include boundary terms:
rhs(:,1) = rhs(:,1) - u_solution(Xindices, 1)/h^2;
rhs(:,m) = rhs(:,m) - u_solution(Xindices,m+2)/h^2;
rhs(1,:) = rhs(1,:) - u_solution(1,Yindices)/h^2;
rhs(m,:) = rhs(m,:) - u_solution(m+2,Yindices)/h^2;

form_rhs = @(m, f, u) 0;

rhs = reshape(rhs, m*m, 1);

%%
u_est = A\rhs;


%%
[X,Y]=meshgrid(0:h:1,0:h:1);
figure;surf(X,Y,u_exact(X,Y));
%%
range = 0:h:1;
[X,Y]=meshgrid(range, range);
figure;surf(X,Y,u_solution);
%%
range = 0:h:1;
range = range(2:end-1);
[X,Y]=meshgrid(range, range);
figure;surf(X,Y,reshape(u_est, m, m));