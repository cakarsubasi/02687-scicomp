a = 0.5;
u_param = @(eps) @(t, x) -tanh((x + 0.5 - t)/(2*eps)) + 1;
u_actual = u_param(a);
%%
m = 50;
points = linspace(-1, 1, m);
U_initial = u_actual(0.0, points);
##for t=0:0.05:2.5
##    U_actual = u_actual(t, points);
##    plot(points, U_actual);
##    pause(0.1);
##end
%%

%%
G = U_initial;
G(2:end) = 0;
maxstep = 500;
plotinterval = 1;
to_plot = true;

err = 0;
err_m = zeros([1 maxstep]);


time = 0.0;
U_new = U_initial;

for i = 1:maxstep
    [U_new, k] = advection_solver(U_new, G, m, a);
    time = time + k;
    U_real = u_actual(time, points);
    % enforce BC
    G(2:end-1) = 0;
    G(1) = U_real(1);
    G(end) = U_real(end);
    err = step_err(U_real, U_new);
    err_m(i) = err;

    if mod(i, plotinterval) == 0 && to_plot
        plot(points, U_new)
        hold on;
        plot(points, U_real)
        legend(["estimate", "real"]);
        title(sprintf("time: %f, err: %f", time, err));
        hold off;
        pause(0.001);
    end
end

#hold on;
#plot(err_m);

