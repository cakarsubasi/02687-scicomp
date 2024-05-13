a = 0.5;
u_param = @(eps) @(t, x) -tanh((x + 0.5 - t)/(2*eps)) + 1;
u_actual = u_param(a);

start = 3;
total = 8;

mvals = 2.^(linspace(start, start+total-1, total))+1;
errors = zeros(total, 1);

for iter = 1:total
m = mvals(iter);

points = linspace(-1, 1, m);
U_initial = u_actual(0.0, points);
figure;

G = U_initial;
G(1:end) = 0;

maxstep = 200000;
plotinterval = 1000;
to_plot = true;
wait = 0;

tmax = 1.0;

time = 0.0;
U_new = U_initial;
for i = 1:maxstep
    maxk = tmax - time;
    if maxk == 0
        plotmain(points, U_new, U_real, err, i, time, wait)
        break
    end
    [U_new, k] = advection_solver(U_new, G, m, a, maxk);
    time = time + k;
    U_real = u_actual(time, points);
    % enforce BC
    G(2:end-1) = 0;
    G(1) = U_real(1);
    G(end) = U_real(end);

    err = step_err(U_real, U_new, m);

    if mod(i, plotinterval) == 0 && to_plot
        plotmain(points, U_new, U_real, err, i, time, wait)
    end
end

errors(iter) = err;
end
%%
figure;
loglog(1./mvals, errors, "o-");
title("Convergence test: advection-diffusion")
xlabel("h");
ylabel("error (Mean Absolute Error)");
grid on;
fig = gcf;
exportgraphics(fig, "ex4_1_convergence.pdf");

%%
function plotmain(points, U_new, U_real, err, i, time, wait)
        plot(points, U_new)
        hold on;
        xlabel("x");
        ylabel("t");
        plot(points, U_real)
        legend(["estimate", "real"]);
        title(sprintf("iter: %d, time: %f, err: %f", i, time, err));
        hold off;
        pause(wait);
end

function err = step_err(est, act, n)
    act = reshape(act, [n, 1]);
    est = reshape(est, [n, 1]);
    err = sum(abs(act - est), "all") / n;
end