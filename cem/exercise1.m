delta = 0.0002;
rhs = @(t, y) y.^2 - y.^3;
%rhs_u = @(t, y) 2*y - 3*y^2;
reltol = 1.e-6;
abstol = 1.e-4;

[T, Y] = ode_solver(rhs, [0 2/delta], delta, reltol, abstol);
plot(T, Y, "o-");
hold on;
%%
opts = odeset('RelTol', reltol, "AbsTol", abstol);
[Tref, Yref] = ode23s(rhs, [0 2/delta], delta, opts);
plot(Tref, Yref, "o-");

[Tref2, Yref2] = ode45(rhs, [0 2/delta], delta, opts);
%plot(Tref2, Yref2, "o-");

legend(["ours", "ode23s", "ode45"])