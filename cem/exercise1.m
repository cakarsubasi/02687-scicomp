delta = 0.0002;
rhs = @(t, y) y^2 - y^3;
%rhs_u = @(t, y) 2*y - 3*y^2;
%opts = odeset('RelTol', 1.e-6);
%ode23s(rhs, [0 2/delta], delta, opts);
[T, Y] = ode_solver(rhs, [0 2/delta], delta, 1.e-6, 1.e-3);
plot(T, Y, "o-");