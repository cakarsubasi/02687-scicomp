% exact solution
u_exact = @(x, y) cos(4*pi*x.*y) + sin(4*pi*(x + y));
% u_xx + u_yy (right hand side function)
f = @(x, y) - 32*pi^2*sin(4*pi*(x + y)) - 16*x.^2*pi^2.*cos(4*pi*x.*y) - 16*y.^2*pi^2.*cos(4*pi*x.*y);
del_f = @(x, y) 1024*pi^4*sin(4*pi*(x + y)) - 64*pi^2*cos(4*pi*x.*y) + 256*x.^4.*pi^4.*cos(4*pi*x.*y) + 256*y.^4.*pi^4.*cos(4*pi*x.*y) + 512*x.*y.*pi^3.*sin(4*pi*x.*y) + 512*x.^2.*y.^2.*pi^4.*cos(4*pi*x.*y);
%%
m = 50;
bc = [0 1 0 1];
[rhs, u_solution, h] = makerhs(f, "9-point", bc, m, u_exact);
A = poisson9(m);
u_est = reshape(A\rhs, m, m);

range = 0:h:1;
range = range(2:end-1);
[X,Y]=meshgrid(range, range);
subplot(1,3,2);
surf(X,Y,reshape(u_est, m, m));
%% convergence test
% 9-point
err = zeros(100, 1);
errX = zeros(100, 1);
for m = 3:103
    bc = [0 1 0 1];
    [rhs, u_solution, h] = makerhs(f, "9-point", bc, m, u_exact);
    A = poisson9(m);
    u_est = reshape(A\rhs, m, m);

    err(m-2) = sqerr(u_est, u_solution);
    errX(m-2) = h;
end
loglog(errX,err);
hold on;
%%
% 5-point
err = zeros(100, 1);
errX = zeros(100, 1);
for m = 3:103
    bc = [0 1 0 1];
    [rhs, u_solution, h] = makerhs(f, "5-point", bc, m, u_exact);
    A = poisson5(m);
    u_est = reshape(A\rhs, m, m);

    err(m-2) = sqerr(u_est, u_solution);
    errX(m-2) = h;
end
loglog(errX,err);
%%
% Amult
err = zeros(100, 1);
errX = zeros(100, 1);
for m = 3:103
    bc = [0 1 0 1];
    [rhs, u_solution, h] = makerhs(f, "5-point", bc, m, u_exact);
    u_est = reshape(Amult(rhs, m), m, m);

    err(m-2) = sqerr(u_est, u_solution);
    errX(m-2) = h;
end
loglog(errX,err);
legend("9-point", "5-point", "amult")

%% Create figure
f = figure;
f.Name = "9-Point Poisson";
f.Position = [10 450 1800 450];

%
range = 0:h:1;
range = range(2:end-1);
[X,Y]=meshgrid(range, range);

subplot(1,3,1);
surf(X,Y,u_solution(2:m+1, 2:m+1));
%
range = 0:h:1;
range = range(2:end-1);
[X,Y]=meshgrid(range, range);
subplot(1,3,2);
surf(X,Y,reshape(u_est, m, m));
% Errors
err = sqrt((u_solution(2:m+1, 2:m+1) - reshape(u_est, m, m)).^2);
range = 0:h:1;
range = range(2:end-1);
[X,Y]=meshgrid(range, range);
subplot(1,3,3);
surf(X,Y,reshape(err, m, m));


%%

% option 1
rhs_err1 = (h^2 / 12) * del_f(Xint, Yint);

% built in operator doesn't make the mistake
% del2 estimates 0.25*(u_xx+y_yy), so we need to multiply by 4
rhs_err2 = (4*h^2/12)*del2(rhs, h);


%%
function err = sqerr(u_est, u_sol)
    err = sum((u_sol(2:end-1, 2:end-1) - u_est).^2, "all") / numel(u_est);
end

function err = abserr(u_est, u_sol)
    err = sum(abs(u_sol(2:end-1, 2:end-1) - u_est), "all") / numel(u_est);
end