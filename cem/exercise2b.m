% exact solution
u_exact = @(x, y) cos(4*pi*x.*y) + sin(4*pi*(x + y));
% u_xx + u_yy (right hand side function)
f = @(x, y) - 32*pi^2*sin(4*pi*(x + y)) - 16*x.^2*pi^2.*cos(4*pi*x.*y) ...
    - 16*y.^2*pi^2.*cos(4*pi*x.*y);
% analytically calculated del_f
del_f = @(x, y) 1024*pi^4*sin(4*pi*(x + y)) - 64*pi^2*cos(4*pi*x.*y) ...
    + 256*x.^4.*pi^4.*cos(4*pi*x.*y) + 256*y.^4.*pi^4.*cos(4*pi*x.*y) ...
    + 512*x.*y.*pi^3.*sin(4*pi*x.*y) + 512*x.^2.*y.^2.*pi^4.*cos(4*pi*x.*y);
%%
m = 50;
bc = [0 1 0 1];

[rhs, u_solution, h] = makerhs(f, "9-point", bc, m, u_exact);
A = poisson9(m);
u_est = reshape(A\rhs, m, m);

range = 0:h:1;
range = range(2:end-1);
[X,Y]=meshgrid(range, range);
figure();
surf(X,Y,reshape(u_est, m, m));
title("Actual solution");
xlabel("x");
ylabel("y");

m_vals = 2.^(1:8); %3:10:103;
err = zeros(size(m_vals));
errX = zeros(size(m_vals));
%% Graphs
m = 50;
% 5-point
[rhs, u_solution, h] = makerhs(f, "5-point", bc, m, u_exact);
A = poisson5(m);
u_est = reshape(A\rhs, m, m);
plot_estimate_and_error("5-point", u_est, u_solution, m, h);

% 9-point del2
[rhs, u_solution, h] = makerhs(f, "9-point", bc, m, u_exact);
A = poisson9(m);
u_est = reshape(A\rhs, m, m);
plot_estimate_and_error("9-point del2", u_est, u_solution, m, h);

% 9-point analytic
[rhs, u_solution, h] = makerhs2(f, bc, m, u_exact, del_f);
A = poisson9(m);
u_est = reshape(A\rhs, m, m);
plot_estimate_and_error("9-point analytic", u_est, u_solution, m, h);

% PCG
[rhs, u_solution, h] = makerhs(f, "5-point", bc, m, u_exact);
Afun = poisson5_f(m);
u_est = reshape(-pcg(Afun, rhs, 1e-6, 200), m, m);
plot_estimate_and_error("PCG", u_est, u_solution, m, h);
%%
omega = 0.4; %2/3;
iters = 100;
% Jacobi
[F, u_solution, ~] = makerhs(f, "5-point", bc, m, u_exact);
U = zeros(m*m, 1);
for i = 1:iters
    U = smooth(U,omega,m,F);
end
u_est = reshape(U, m, m);
plot_estimate_and_error("Jacobi: \omega = 0.4", u_est, u_solution, m, h);
omega = 0.66; %2/3;
% Jacobi
[F, u_solution, ~] = makerhs(f, "5-point", bc, m, u_exact);
U = zeros(m*m, 1);
for i = 1:iters
    U = smooth(U,omega,m,F);
end
u_est = reshape(U, m, m);
plot_estimate_and_error("Jacobi: \omega = 0.66", u_est, u_solution, m, h);
omega = 0.99; %2/3;
% Jacobi
[F, u_solution, ~] = makerhs(f, "5-point", bc, m, u_exact);
U = zeros(m*m, 1);
for i = 1:iters
    U = smooth2(U,omega,m,F);
end
u_est = reshape(U, m, m);
plot_estimate_and_error("Jacobi: \omega = 0.80", u_est, u_solution, m, h);


%% convergence test
% 9-point
i = 1;
for m = m_vals
    bc = [0 1 0 1];
    [rhs, u_solution, h] = makerhs(f, "9-point", bc, m, u_exact);
    A = poisson9(m);
    u_est = reshape(A\rhs, m, m);

    err(i) = abserr(u_est, u_solution);
    errX(i) = h;
    i = i + 1;
end
plot_convergence("9-point direct - del2", errX, err)
%%
% 9-point analytic
i = 1;
for m = m_vals
    bc = [0 1 0 1];
    [rhs, u_solution, h] = makerhs2(f, bc, m, u_exact, del_f);
    A = poisson9(m);
    u_est = reshape(A\rhs, m, m);

    err(i) = abserr(u_est, u_solution);
    errX(i) = h;
    i = i + 1;
end
plot_convergence("9-point direct - analytic", errX, err);
%%
% 5-point
i = 1;
for m = m_vals
    [rhs, u_solution, h] = makerhs(f, "5-point", bc, m, u_exact);
    A = poisson5(m);
    u_est = reshape(A\rhs, m, m);

    err(i) = abserr(u_est, u_solution);
    errX(i) = h;
    i = i + 1;
end
plot_convergence("5-point direct", errX, err);
%%
% pcg
i = 1;
for m = m_vals
    bc = [0 1 0 1];
    [rhs, u_solution, h] = makerhs(f, "5-point", bc, m, u_exact);
    Afun = poisson5_f(m);
    u_est = reshape(-pcg(Afun, rhs, 1e-6, 20000), m, m);

    err(i) = abserr(u_est, u_solution);
    errX(i) = h;
    i = i + 1;
end
plot_convergence("pcg", errX, err);

%%
% Jacobi smoothing
i = 1;
maxiter = 20000;
for m = m_vals
    bc = [0 1 0 1];
    [F, u_solution, h] = makerhs(f, "5-point", bc, m, u_exact);
    tol = 1e-6;
    U = zeros(m*m, 1);
    for j = 1:maxiter
        U = smooth(U,omega,m,F);
        R = F + Amult(U, m);
        res = sum(abs(R), "all");
        if res < tol
            fprintf("Jacobi converged in %d iterations, residual: %f\n", j, res);
            break;
        elseif j == maxiter
            fprintf("Jacobi did not converge in %d iterations, residual: %f\n", j, res);
        end
    end

    u_est = reshape(U, m, m);

    err(i) = abserr(u_est, u_solution);
    errX(i) = h;
    i = i + 1;
end
plot_convergence("Jacobi", errX, err);

%% Create figure
fig = figure;
fig.Name = "9-Point Poisson";
fig.Position = [10 450 1800 450];

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
function err = sqerr(u_est, u_sol)
    err = sum((u_sol(2:end-1, 2:end-1) - u_est).^2, "all") / numel(u_est);
end

function err = abserr(u_est, u_sol)
    err = sum(abs(u_sol(2:end-1, 2:end-1) - u_est), "all") / numel(u_est);
end
%% plotting routines
function plot_estimate_and_error(name, u_est, u_solution, m, h)

f = figure;
f.Name = name;
f.Position = [10 450 1200 380];

%
range = 0:h:1;
range = range(2:end-1);
[X,Y]=meshgrid(range, range);
subplot(1,2,1);
surf(X,Y,reshape(u_est, m, m));
title("Estimated solution");
xlabel("x");
ylabel("y");
% Errors
err = sqrt((u_solution(2:m+1, 2:m+1) - reshape(u_est, m, m)).^2);
range = 0:h:1;
range = range(2:end-1);
[X,Y]=meshgrid(range, range);
subplot(1,2,2);
surf(X,Y,reshape(err, m, m));
title("Error");
xlabel("x");
ylabel("y");
end

function plot_convergence(name, errX, err)
f = figure;
f.Name = name;
loglog(errX,err);
grid on;
title(sprintf("Convergence graph: %s", name));
xlabel("h");
ylabel("error");
end
