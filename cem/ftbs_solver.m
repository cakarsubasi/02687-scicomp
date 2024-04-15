function [U_new, k] = ftbs_solver(U, G, m, a)
%FTBS_SOLVER our ftbs solver
%Solves problems of the form u_t(x,t) + a*u_x(x, t) = 0
%   * U: [m, 1] previous iteration
%   * G: [m, 1] boundary conditions, only the first or last element is used
%   * m: size of vector
%   * a: parameter a (determines which sided variant will be used)
%
%   * U_new: new iteration
%   * k: step taken

h = 2/(m+1);

% pick k values
% left handed variant
if a > 0
    k = h^2 / a;
    cons = k*(-a)/h;
% right handed variant
else
    k = h^2 / -a;
    cons = -k*a/h;
end

U = reshape(U, [m, 1]);
G = reshape(G, [m, 1]);
% left handed variant
if a > 0
    e = ones(m, 1);
    z = zeros(m, 1);
    A = spdiags([-1*e 1*e z], [-1 0 1], m, m)*cons;
    G(2:end) = 0;
% right handed variant
else
    e = ones(m, 1);
    z = zeros(m, 1);
    A = spdiags([z -1*e 1*e], [-1 0 1], m, m)*cons;
    G(1:end-1) = 0;
    G(end) = -G(end);
end

U_new = U + A*U - cons*G;
U_new = reshape(U_new, [m, 1]);
end

