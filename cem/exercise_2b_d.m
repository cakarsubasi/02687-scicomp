% exact solution
u_exact = @(x, y) cos(4*pi*x.*y) + sin(4*pi*(x + y));
% u_xx + u_yy (right hand side function)
f = @(x, y) - 32*pi^2*sin(4*pi*(x + y)) - 16*x.^2*pi^2.*cos(4*pi*x.*y) - 16*y.^2*pi^2.*cos(4*pi*x.*y);
del_f = @(x, y) 1024*pi^4*sin(4*pi*(x + y)) - 64*pi^2*cos(4*pi*x.*y) + 256*x.^4.*pi^4.*cos(4*pi*x.*y) + 256*y.^4.*pi^4.*cos(4*pi*x.*y) + 512*x.*y.*pi^3.*sin(4*pi*x.*y) + 512*x.^2.*y.^2.*pi^4.*cos(4*pi*x.*y);
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
Xint = X(Xindices,Yindices);       % interior points
Yint = Y(Xindices,Yindices);

A = poisson9(m);

% solution:
u_solution = u_exact(X, Y);
%%
rhs = f(Xint, Yint);
rhs2 = rhs + (h^2 / 12) * del_f(Xint, Yint);

rhs_all = f(X, Y);
rhs_err = 0.5 * (h^2 / 12) * ... % why the 0.5?
    (-4 * rhs_all(Xindices, Yindices) ...
    + rhs_all(Xindices+1, Yindices+1) ...
    + rhs_all(Xindices+1, Yindices-1) ...
    + rhs_all(Xindices-1, Yindices+1) ...
    + rhs_all(Xindices-1, Yindices-1))/h^2;

%figure;
%surf(0.5 * rhs_err);
%figure;
%surf((h^2 / 12) * del_f(Xint, Yint));
%surf(del_f(Xint, Yint) - 2*rhs_err)

%%

rhs = rhs_all(Xindices, Yindices) + rhs_err;


% adjust the rhs to include boundary terms:
rhs(:,1) = rhs(:,1) - (4*u_solution(Xindices, 1) + u_solution(Xindices-1, 1) + u_solution(Xindices+1, 1))/6/h^2;  % top
rhs(:,m) = rhs(:,m) - (4*u_solution(Xindices,m+2) + u_solution(Xindices-1,m+2) + u_solution(Xindices+1,m+2))/6/h^2; % bottom
rhs(1,:) = rhs(1,:) - (4*u_solution(1,Yindices) + u_solution(1,Yindices-1) + u_solution(1,Yindices+1))/6/h^2;   % left
rhs(m,:) = rhs(m,:) - (4*u_solution(m+2,Yindices) + u_solution(m+2,Yindices-1) + u_solution(m+2,Yindices+1))/6/h^2; % right
% fix the corners
rhs(1,1)      = rhs(1, 1)     + u_solution(1, 1)    /6/h^2;
rhs(1,end)    = rhs(1, end)   + u_solution(1, end)  /6/h^2;
rhs(end,1)    = rhs(end, 1)   + u_solution(end, 1)  /6/h^2;
rhs(end, end) = rhs(end, end) + u_solution(end, end)/6/h^2;
%%
rhs = reshape(rhs, m*m, 1);

u_est = A\rhs;

%%
%range = 0:h:1;
%[X,Y]=meshgrid(range, range);
%figure;surf(X,Y,u_solution);
%% Create figure
f = figure;
f.Name = "9-Point Poisson";
f.Position = [10 450 1800 450];

%%
range = 0:h:1;
range = range(2:end-1);
[X,Y]=meshgrid(range, range);

subplot(1,3,1);
surf(X,Y,u_solution(2:m+1, 2:m+1));
%%
range = 0:h:1;
range = range(2:end-1);
[X,Y]=meshgrid(range, range);
subplot(1,3,2);
surf(X,Y,reshape(u_est, m, m));
%% Errors
err = sqrt((u_solution(2:m+1, 2:m+1) - reshape(u_est, m, m)).^2);
range = 0:h:1;
range = range(2:end-1);
[X,Y]=meshgrid(range, range);
subplot(1,3,3);
surf(X,Y,reshape(err, m, m));