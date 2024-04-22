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

err = 0;
err_m = zeros([1 maxstep]);

maxstep = 200;
plotinterval = 1;
to_plot = true;

time = 0;
U_new = U_initial;
for i = 1:maxstep
    [U_new, k] = ftcs_solver(U_new, G, m, eps, 1.04);
    time = time + k;
    U_real = u_actual(time, points);
    G = [u_actual(time, points(1)), u_actual(time, points(end))];

    err = step_err(U_real, U_new);
    err_m(i) = err;

    if mod(i, plotinterval) == 0 && to_plot
        plot([G(1); U_new; G(2)]);
        hold on;
        plot(U_real)
        legend(["estimate", "real"]);
        title(sprintf("time: %f, err: %f", time, err));
        hold off;
        pause(0.1);
    end
end

function err = step_err(est, act)
    err = sum(abs(act - est), "all") / length(est);
end