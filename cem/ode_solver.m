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

% apply the RK method

A1 = [0  0 0;
     0.5 0 0;
     0.5 0.5 0.0];
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
    
    if norm(err) > reps*norm(y_new2) + aeps
        step = step/2;
    else
        t_n = t_n + step;
        step = tend - t_n;
        y_n = y_new1;
        YOUT = cat(2, YOUT, y_new1);
        TOUT = cat(2, TOUT, t_n);
    end

end


end

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
