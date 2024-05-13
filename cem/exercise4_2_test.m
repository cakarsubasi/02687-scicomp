exact = -152.00516;

start = 3;
total = 10;

mvals = 2.^(linspace(start, start+total-1, total))+1;
estimates = zeros(total, 1);

for iter = 1:total
m = mvals(iter);

points = linspace(-1, 1, m);
U_initial = -sin(pi*points);
a = 0.01/pi;
figure;
plot(U_initial);

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

plotmain(U_new, U_initial, i, time, wait)

figure;
U_x = (- U_new(1:end-2) + U_new(3:end))/(4/(m+1));
wavebreak = min(U_x);
plot(points(2:end-1), U_x);
title(sprintf("$$U_x$$ at 0: %f", wavebreak), "Interpreter", "latex");

disp(wavebreak);

estimates(iter) = wavebreak;
end
%%

figure;
loglog(1./mvals, (exact - estimates).^2, "o-");
title("Convergence test: breaking wave")
xlabel("h");
ylabel("error (MSE)");
grid on;
fig = gcf;
exportgraphics(fig, "ex4_convergence2.pdf");

%%
figure;
loglog(1./mvals(2:end), diff(estimates).^2, "o-");
title("Residuals: breaking wave")
xlabel("h");
ylabel("residuals squared");
grid on;
fig = gcf;
exportgraphics(fig, "ex4_residuals2.pdf");


function plotmain(U_new, U_initial, i, time, wait)
        plot(U_new)
        hold on;
        plot(U_initial)
        legend(["estimate", "initial"]);
        title(sprintf("iter: %d, time: %f", i, time));
        hold off;
        pause(wait);
end