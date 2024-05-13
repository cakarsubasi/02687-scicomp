a = 0.5;
u_param = @(eps) @(t, x) -tanh((x + 0.5 - t)/(2*eps)) + 1;
u_actual = u_param(a);
%%
m = 100;
points = linspace(-1, 1, m);
U_initial = u_actual(0.0, points);
for t=0:0.1:0
    U_actual = u_actual(t, points);
    plot(U_actual);
    pause(0.1);
end

%%
G = U_initial;
G(2:end) = 0;

maxstep = 200000;
plotinterval = 50;
to_plot = true;
wait = 0.01;

err = 0;
err_m = zeros([1 maxstep]);
tmax = 1.0;

time = 0.0;
U_new = U_initial;
for i = 1:maxstep
    maxk = tmax - time;
    if maxk == 0
        plotmain(points, U_new, U_real, err, i, time, wait);
        break
    end
    [U_new, k] = advection_solver(U_new, G, m, a, maxk);
    time = time + k;
    U_real = u_actual(time, points);
    % enforce BC
    G(2:end-1) = 0;
    G(1) = U_real(1);
    G(end) = U_real(end);

    err = step_err(U_real, U_new);
    err_m(i) = err;

    if mod(i, plotinterval) == 0 && to_plot
        plotmain(points, U_new, U_real, err, i, time, wait);
    end
end

%hold on;
%plot(err_m);
function plotmain(points, U_new, U_initial, err, i, time, wait)
        plot(points, U_new)
        hold on;
        xlabel("x");
        ylabel("t");
        plot(points, U_initial)
        legend(["estimate", "real"]);
        title(sprintf("iter: %d, time: %f, err: %f", i, time, err));
        hold off;
        pause(wait);
end


function err = step_err(est, act)
    err = sum(abs(act - est), "all") / length(est);
end