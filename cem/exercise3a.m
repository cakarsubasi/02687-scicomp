%%
% exact solution
u_exact = @(x, y) cos(4*pi*x.*y) + sin(4*pi*(x + y));
% right hand side function
f = @(x, y) - 32*pi^2*sin(4*pi*(x + y)) - 16*x.^2*pi^2.*cos(4*pi*x.*y) - 16*y.^2*pi^2.*cos(4*pi*x.*y);
%% BCs
a = 0;
b = 1;
c = 0;
d = 1;

m = 50;
h = (b-a)/(m+1);

X = linspace(a, b, m+2); % all points
Y = linspace(c, d, m+2);
[X,Y] = meshgrid(X,Y);
Xindices = 2:m+1; % interior points
Yindices = 2:m+1;
Xint = X(Xindices,Yindices); % interior points
Yint = Y(Xindices,Yindices);

% solution:
u_solution = u_exact(X, Y);
%%
rhs = f(Xint, Yint);

% adjust the rhs to include boundary terms:
rhs(:,1) = rhs(:,1) - u_solution(Xindices, 1)/h^2;
rhs(:,m) = rhs(:,m) - u_solution(Xindices,m+2)/h^2;
rhs(1,:) = rhs(1,:) - u_solution(1,Yindices)/h^2;
rhs(m,:) = rhs(m,:) - u_solution(m+2,Yindices)/h^2;

%%
%rhs = reshape(rhs, m*m, 1);
%u_est = Amult(rhs, m);

% 3.2
%rhs = f(X, Y);
u_curr = zeros(m, m);
for i = 1:2
u_curr = jacobis(u_curr, 0.67, m, rhs);
end
u_est = u_curr;

%%
range = 0:h:1;
[X,Y]=meshgrid(range, range);
figure;surf(X,Y,u_solution);
%%
range = 0:h:1;
range = range(2:end-1);
[X,Y]=meshgrid(range, range);
figure;surf(X,Y,u_solution(2:m+1, 2:m+1));
%%
range = 0:h:1;
range = range(2:end-1);
[X,Y]=meshgrid(range, range);
figure;surf(X,Y,reshape(u_est, m, m));
%% Errors
err = sqrt((u_solution(2:m+1, 2:m+1) - reshape(u_est, m, m)).^2);
range = 0:h:1;
range = range(2:end-1);
[X,Y]=meshgrid(range, range);
figure;surf(X,Y,reshape(err, m, m));