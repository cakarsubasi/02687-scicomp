%%
% This file just contains a differentiation to validate paper
% differentiation
syms x y
u = sin(4*pi*(x+y)) + cos(4*pi*x.*y);
u_xx = diff(u, x, 2);
u_yy = diff(u, y, 2);
f_ex = u_xx + u_yy;

del_f_ex = diff(f_ex, x, 2) + diff(f_ex, y, 2);

u_exact = @(a, b) subs(u, [x y], [a b]);
f = @(a, b) subs(f_ex, [x y], [a b]);
del_f = @(a, b) subs(del_f_ex, [x y], [a, b]); 