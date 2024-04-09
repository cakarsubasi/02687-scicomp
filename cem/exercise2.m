alpha = [1; 4; 16];
a = 1;
b = 0;
eps = 0.1;

u_actual = @(t, x) sum(exp(-eps*alpha.^2.*t).*(a.*cos(alpha.*x)+b.*sin(alpha.*x)), 1);

%%
m = 64;
points = linspace(-1, 1, m);
U_initial = u_actual(0.0, points);
G = U_initial;
G(2:end-1) = 0;

h = 2/(m+1);
k = h^2 / (2*eps);
% pick h and k values
cons = eps*k/h^2;

e = ones(m, 1);
A = spdiags([1*e -2*e 1*e], [-1 0 1], m, m)*cons;

U_initial = U_initial';
time = 0.0;
for i = 1:200
    U_new = U_initial + A*U_initial + cons*G';
    time = time + k;
    plot(U_new)
    G = u_actual(time, points);
    G(2:end-1) = 0;
    U_initial = U_new;
    pause(0.5);
end
%[U_new, k] = ftcs_solver(U_initial, G, m, eps);