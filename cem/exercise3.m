f = @(y) sin(2*pi*y);
a = 0.5;
u_param = @(a) @(t, x) f(x - a*t);
u_actual = u_param(a);
%%
m = 200;
points = linspace(-1, 1, m);
U_initial = u_actual(0, points);
G = U_initial;
G(2:end) = 0;

maxstep = 200000;
plotinterval = 2000;
to_plot = true;
wait = 0.01;
tmax = 4/a;

err = 0;
err_m = zeros([1 maxstep]);


time = 0.0;
U_new = U_initial;
figure;
for i = 1:maxstep
    kmax = tmax - time;
    if kmax == 0
        plotmain(points, U_new, U_real, G, err, i, time, wait);
        break
    end
    [U_new, k] = ftbs_solver(U_new, G, m, a, kmax);

    time = time + k;
    U_real = u_actual(time, points);
    % enforce BC
    G(2:end) = 0;
    G(1) = U_new(end);

    err = step_err(U_real, U_new, m);
    err_m(i) = err;

    if mod(i, plotinterval) == 0 && to_plot
        plotmain(points, U_new, U_real, G, err, i, time, wait);
    end
end

exportgraphics(gcf, "ex3_t4_200.pdf");

%%
function plotmain(points, U_new, U_real, G, err, iter, time, wait)
        plot(points, U_new);
        hold on;
        plot(points, U_real)
        legend(["estimate", "real"]);
        title(sprintf("time: %f, err: %f, iter: %d", time, err, iter));
        xlabel("x");
        ylabel("U")
        hold off;
        pause(wait);
end

function err = step_err(est, act, n)
    act = reshape(act, [n, 1]);
    est = reshape(est, [n, 1]);
    err = sum(abs(act - est), "all") / n;
end