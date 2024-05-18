% Exercise 4 - Nonlinear Advection Equation solver on the shockwave
m = 4097;
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
plotmain(points, U_new, U_initial, i, time, wait)
fig = gcf;
exportgraphics(fig, "ex4_2_example.pdf");

figure;
%U_x = diff(U_new)/(2/(m+1));
U_x = (- U_new(1:end-2) + U_new(3:end))/(4/(m+1));
wavebreak = min(U_x);
plot(points(2:end-1), U_x);
title(sprintf("$$\\partial_x U; \\quad \\partial_x U |_{x=0} = %f$$", wavebreak), "Interpreter", "latex");
ylabel(sprintf("$$\\partial_x U$$"), "Interpreter", "latex")
xlabel(sprintf("$$x$$"), "Interpreter", "latex")
fig = gcf;
exportgraphics(fig, "ex4_2_shockwave.pdf");
disp(wavebreak);

function plotmain(points, U_new, U_initial, i, time, wait)
        plot(points, U_new)
        hold on;
        plot(points, U_initial)
        legend(["estimate", "initial"]);
        title(sprintf("iter: %d, time: %f", i, time));
        xlabel("x");
        ylabel(sprintf("U(x,%f)", time));
        hold off;
        pause(wait);
end