%%%%%%%%%%%%%%%%% Question 2%%%%%%%%%%%%%%%

##a = 0.5;
##u_param = @(eps) @(t, x) -tanh((x + 0.5 - t)/(2*eps)) + 1;
##u_actual = u_param(a);
##
##m = 50;
##points = linspace(-1, 1, m+1);
##int_points = points(2:(end-1));
##
##U_initial = u_actual(0.0, int_points);
##maxstep = 500;
##plotinterval = 1;
##err = 0;
##time = 0.0;
##
##U_new = U_initial;
##gL = u_actual(time, -1);
##gR = u_actual(time, 1);
##
##for i = 1:maxstep
##    [U_new, k] = advection_solver(U_new, gL, gR, m, a);
##    time = time + k;
##    U_real = u_actual(time, int_points);
##    % enforce BC
##    gL = u_actual(time, -1);
##    gR = u_actual(time, 1);
##    err = step_err(U_real, U_new);
##
##    if mod(i, plotinterval) == 0
##        plot(int_points, U_new)
##        hold on;
##        plot(int_points, U_real)
##        legend(["estimate", "real"]);
##        title(sprintf("time: %f, err: %f", time, err));
##        hold off;
##        pause(0.001);
##    end
##end

%%%%%%%%%%%%%%%%% Question 3 %%%%%%%%%%%%%%%

##a = 0.01/pi;
##nu_actual = @(t, x) -sin(pi*x);
##
##m = 2000;
##points = linspace(-1, 1, m+1);
##int_points = points(2:(end-1));
##
##U_initial = nu_actual(0.0, int_points);
##maxstep = 5000;
##plotinterval = 1;
##err = 0;
##time = 0.0;
##tf = 1.6037/pi;
##
##U_new = U_initial;
##gL = 0;
##gR = 0;
##
##for i = 1:maxstep
##    [U_new, k] = advection_solver(U_new, gL, gR, m, a);
##    time = time + k;
##    k
##    if time > tf
##      break;
##    endif
##
##    if mod(i, plotinterval) == 0
##        plot(int_points, U_new)
##        title(sprintf("time: %f", time));
##        hold off;
##        pause(0.00001)
##    end
##end
##
##h = 2/(m+1);
##Uend = [0 U_new' 0];
##Uprime = (Uend(1:(end-2))-Uend(3:end))/(2*h);
##Uprime(m/2)

%%%%%%%%%%%%%%%%% Question 4 %%%%%%%%%%%%%%%

a = 0.01/pi;
nu_actual = @(t, x) -sin(pi*x);

m = 2000;
points = linspace(-1, 1, m+1);
int_points = points(2:(end-1));
dist = diff(points);

U_initial = nu_actual(0.0, int_points);
maxstep = 1000;
plotinterval = 1;
err = 0;
time = 0.0;
tf = 1.6037/pi;

U_new = U_initial;
gL = 0;
gR = 0;

for i = 1:maxstep
    [U_new, k] = adapted_mesh_solver(U_new, gL, gR, m, a, dist);
    time = time + k;
    if time > tf
      break;
    endif

    if mod(i, plotinterval) == 0
        plot(int_points, U_new)
        title(sprintf("time: %f", time));
        hold off;
        pause(0.001)
    end
end

U_x = diff(U_new)/(2/(m+1));





