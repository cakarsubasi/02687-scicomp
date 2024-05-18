% plotting scripts for exercise 2

alpha = [1; 4; 16];
a = 1;
b = 0;
eps = 0.1;

u_actual = @(t, x) sum(exp(-eps*alpha.^2.*t).*(a.*cos(alpha.*x)+b.*sin(alpha.*x)), 1);

%%
m = 64;
points = linspace(-1, 1, m+2);
U_initial = u_actual(0.0, points(2:end-1));
G = [u_actual(0, points(1)), u_actual(0, points(end))];

% set this above 1.0 to make the method unstable
step_fac = 1.0;

err = 0;
err_m = zeros([1 maxstep]);

maxstep = 600;
plotinterval = 2;
to_plot = true;
tmax = 1;
wait = 0.5;

time = 0;
U_new = U_initial;
for i = 1:maxstep
    kmax = tmax - time;
    if kmax == 0
        break
    end
    [U_new, k] = ftcs_solver(U_new, G, m, eps, step_fac, kmax);
    time = time + k;
    U_real = u_actual(time, points);
    G = [u_actual(time, points(1)), u_actual(time, points(end))];

    err = step_err(U_real, U_new);
    err_m(i) = err;

    if mod(i, plotinterval) == 0 && to_plot
        plotmain(points, U_new, U_real, G, err, i, time, wait);
    end
end

plotmain(points, U_new, U_real, G, err, i, time, wait);

function plotmain(points, U_new, U_real, G, err, iter, time, wait)
        plot(points, [G(1); U_new; G(2)]);
        hold on;
        plot(points, U_real)
        legend(["estimate", "real"]);
        title(sprintf("time: %f, err: %f, iter: %d", time, err, iter));
        xlabel("x");
        ylabel("U")
        hold off;
        pause(wait);
end

function err = step_err(est, act)
    err = sum(abs(act - est), "all") / length(est);
end