f = @(y) sin(2*pi*y);
a = 0.5;
u_param = @(a) @(t, x) f(x - a*t);
u_actual = u_param(a);
%%
start = 3;
total = 8;

mvals = 2.^(linspace(start, start+total-1, total));
errors = zeros(total, 1);
powers = zeros(total, 1);
powers_real = zeros(total, 1);
phases = zeros(total, 1);

for iter = 1:total
m = mvals(iter);

points = linspace(-1, 1, m);
U_initial = u_actual(0, points);
G = U_initial;
G(2:end) = 0;

maxstep = 200000;
plotinterval = 2000;
to_plot = true;
wait = 0.01;
tmax = 1/a;

err = 0;
err_m = zeros([1 maxstep]);


time = 0.0;
U_new = U_initial;
figure;
for i = 1:maxstep
    kmax = tmax - time;
    if kmax == 0
        plotmain(points, U_new, U_real, err, i, time, wait);
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
        plotmain(points, U_new, U_real, err, i, time, wait);
    end
end

errors(iter) = err;
powers(iter) = signal_power(U_new);
powers_real(iter) = signal_power(U_real);
phases(iter) = phase_diff(U_new, U_real);
end
%%
figure;
loglog(1./mvals, errors, "o-");
title("Convergence test: unsteady advection")
xlabel("h");
ylabel("error (Mean Absolute Error)");
grid on;
exportgraphics(gcf, "ex3_convergence.pdf");

figure;
loglog(1./mvals, powers, "o-");
hold on;
loglog(1./mvals, powers_real, "o-");
title("Powers: unsteady advection")
legend(["estimate", "real"]);
xlabel("h");
ylabel("Power");
grid on;
exportgraphics(gcf, "ex3_power.pdf");

figure;
semilogx(1./mvals, powers, "o-");
hold on;
semilogx(1./mvals, powers_real, "o-");
title("Powers: unsteady advection")
legend(["estimate", "real"]);
xlabel("h");
ylabel("Power");
grid on;
exportgraphics(gcf, "ex3_power_linear.pdf");

figure;
semilogx(1./mvals, phases, "o-");
hold on;
title("Phases: unsteady advection")
xlabel("h");
ylabel("Phase-diff");
grid on;
exportgraphics(gcf, "ex3_phase.pdf");


%%
function plotmain(points, U_new, U_real, err, iter, time, wait)
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

function pow = signal_power(signal)
    pow = sum((signal.^2), "all") / length(signal);
end

function dif = phase_diff(signal1, signal2)
   top = dot(signal1, signal2);
   bot = (norm(signal1)*norm(signal2));
   dif = acos(top/bot);
end

function err = step_err(est, act, n)
    act = reshape(act, [n, 1]);
    est = reshape(est, [n, 1]);
    err = sum(abs(act - est), "all") / n;
end