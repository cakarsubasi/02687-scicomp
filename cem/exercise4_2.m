%%
m = 3001;
points = linspace(-1, 1, m);
U_initial = -sin(pi*points);
a = 0.01/pi;
plot(U_initial);
%%
G = U_initial;
G(1:end) = 0;

maxstep = 200000;
plotinterval = 1000;
to_plot = true;
wait = 0;

tmax = 1.6037 / pi;

time = 0.0;
U_new = U_initial;
for i = 1:maxstep
    if time == tmax
        break;
    end
    [U_new, k] = advection_solver2(U_new, G, m, a, time);
    time = time + k;

    if mod(i, plotinterval) == 0 && to_plot
        plotmain(U_new, U_initial, i, time, wait)
    end
end
%%
plotmain(U_new, U_initial, i, time, wait)

figure;
U_x = diff(U_new)/(2/(m+1));
wavebreak = min(U_x);
plot(points(2:end), U_x);
title(sprintf("$$U_x$$ at 0: %f", wavebreak), "Interpreter", "latex");

disp(wavebreak);

function plotmain(U_new, U_initial, i, time, wait)
        plot(U_new)
        hold on;
        plot(U_initial)
        legend(["estimate", "initial"]);
        title(sprintf("iter: %d, time: %f", i, time));
        hold off;
        pause(wait);
end