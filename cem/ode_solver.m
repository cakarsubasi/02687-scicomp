function [TOUT, YOUT] = ode_solver(ODEFUN, TSPAN, Y0, reps, aeps)
%ODE_SOLVER our ode solver
%   * ODEFUN rhs function handle for the ODE in the form f(t, y) which
%     returns a column vector
%   * TSPAN two element vector specifying beginning and end points
%   * Y0 column vector of initial conditions

% apply Picard-Lindeloff to determine step size
tstart = TSPAN(1);
tend = TSPAN(2);

step = tend - tstart;
t_n = tstart;
y_n = Y0;

%
% apply the RK method

A1 = [0  0 0;
     0.5 0 0;
     0.5 0.5 0.0];
B1 = [1/6 2/3 1/6];

A2 = [0.0 0.0 0.0;
      0.5 0.0 0.0;
      -1.0 2.0 0.0];
B2 = B1;
odefun = ODEFUN;

YOUT = [y_n];
TOUT = [tstart];

while t_n < tend

y_new1 = rk23(A1, B1, y_n, t_n, odefun, step);
y_new2 = rk23(A2, B2, y_n, t_n, odefun, step);

err = y_new2 - y_new1;

if norm(err) > reps*norm(y_new2) + aeps
    step = step/2;
else
    t_n = t_n + step;
    y_n = y_new1;
    YOUT = cat(2, YOUT, y_new1);
    TOUT = cat(2, TOUT, t_n);
end

end


end

function y_new = rk23(A, B, y_n, t_n, odefun, step)
%a21 = 0.5;
%a31 = 0.5;
%a32 = 0.5;
a21 = A(2, 1);
a31 = A(3, 1);
a32 = A(3, 2);
%b1 = 0.166666;
%b2 = 0.666666;
%b3 = b1;
b1 = B(1);
b2 = B(2);
b3 = B(3);
%c1 = 0;
%c2 = 0.5;
%c3 = 1;
C = sum(A, 2);
c1 = C(1);
c2 = C(2);
c3 = C(3);


ksi1 = y_n;
ksi2 = y_n + a21*step*odefun(t_n, ksi1);
ksi3 = y_n + a31*step*odefun(t_n, ksi1) + a32*step*odefun(t_n+c2*step, ksi2);

y_new = y_n + step*(b1*odefun(t_n, ksi1) + b2*odefun(t_n+c2*step, ksi2) + b3*odefun(t_n+c3*step, ksi3));
%err = step
end
