f = @(y) sin(2*pi*y);
a = 0.5;
u_param = @(a) @(t, x) f(x - a*t);
u_actual = u_param(a);
%%
m = 64;
points = linspace(-1, 1, m);
U_initial = u_actual(0.0, points);
G = U_initial;
G(2:end) = 0;

maxstep = 2000000;
plotinterval = 50;

time = 0.0;
U_new = U_initial;
for i = 1:maxstep
    [U_new, k] = ftbs_solver(U_new, G, m, a);
    time = time + k;
    U_real = u_actual(time, points);
    G = u_actual(time, points);
    G(2:end-1) = 0;
    if mod(i, plotinterval) == 1
        plot(U_new)
        hold on;
        plot(U_real)
        legend(["estimate", "real"]);
        title(sprintf("%f", time));
        hold off;
        pause(0.01);
    end
end