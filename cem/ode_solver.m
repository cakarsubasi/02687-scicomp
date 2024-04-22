function [TOUT, YOUT] = ode_solver(ODEFUN, TSPAN, Y0, reps, aeps)
%ODE_SOLVER our ode solver
%   * ODEFUN rhs function handle for the ODE in the form f(t, y) which
%     returns a column vector
%   * TSPAN two element vector specifying beginning and end points
%   * Y0 column vector of initial conditions

tstart = TSPAN(1);
tend = TSPAN(2);

step = tend - tstart;
t_n = tstart;
y_n = Y0;

% apply the RK method

A1 = [0    0   0  ;
      0.5  0   0  ;
      0.5  0.5 0.0 ];
B1 = [1/6 2/3 1/6];

A2 = [0.0 0.0 0.0;
      0.5 0.0 0.0;
      -1.0 2.0 0.0];
B2 = B1;

YOUT = y_n;
TOUT = tstart;

while t_n < tend

    y_new1 = rk23(A1, B1, y_n, t_n, ODEFUN, step);
    y_new2 = rk23(A2, B2, y_n, t_n, ODEFUN, step);
    
    err = y_new2 - y_new1;
    
    if norm(err) > abs(norm(y_new2)*reps) + abs(aeps)
        step = step/2;
    else
        t_n = t_n + step;
        step = tend - t_n;
        y_n = y_new1;
        T_star = picard_lindeloef(t_n, tend, y_n, ODEFUN);
        step = T_star;
        YOUT = cat(2, YOUT, y_new1);
        TOUT = cat(2, TOUT, t_n);
    end
end %while

end %function ode_solver

function y_new = rk23(A, B, y_n, t_n, odefun, step)
a21 = A(2, 1);
a31 = A(3, 1);
a32 = A(3, 2);
b1 = B(1);
b2 = B(2);
b3 = B(3);
C = sum(A, 2);
c1 = C(1);
c2 = C(2);
c3 = C(3);


ksi1 = y_n;
ksi2 = y_n + a21*step*odefun(t_n, ksi1);
ksi3 = y_n + a31*step*odefun(t_n, ksi1) + a32*step*odefun(t_n+c2*step, ksi2);

y_new = y_n + step*(b1*odefun(t_n, ksi1) + b2*odefun(t_n+c2*step, ksi2) + b3*odefun(t_n+c3*step, ksi3));
end

function T_star = picard_lindeloef(t0, t1, y0, fun)
T_star = 0;
a = 0.1;
while true

% optimize S
S = max(abs(fun(0, linspace(y0-a, y0+a, 100))));

T_star_new = a/S;
if (T_star_new < T_star) || (t0 + T_star_new >= t1)
    break
else
    T_star = T_star_new;
    a = a + 0.1;
end
end

T_star = min(t1, t0 + T_star);

end
