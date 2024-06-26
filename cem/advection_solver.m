function [U_new, k] = advection_solver(U, G, m, eps, maxk)
% ADVECTION_SOLVER our advection solver
% Solves problems of the form u_t(x,t) + u(x,t)*u_x(x, t) = eps*u_xx(x,t)
%   * U: [m, 1] previous iteration
%   * G: [m, 1] boundary conditions, only the first or last element is used
%   * m: size of vector
%   * eps: parameter eps
%   * MAXK: maximum time step that can be taken
%
%   * U_new: new iteration
%   * k: step taken

h = 2/(m+1);
k = h^2 / (2*eps);
if k > maxk
    k = maxk;
end

e = ones(m, 1);
z = zeros(m, 1);
A1 = spdiags([1*e -2*e 1*e], [-1 0 1], m, m);
A2 = spdiags([-1*e z 1*e], [-1 0 1], m, m);
cons1 = k*eps/h^2;
cons2 = k/(2*h);
U = reshape(U, [m 1]);
G = reshape(G, [m 1]);
G2 = G.*G;
G2(end) = -G2(end);
U_new = U + cons1*A1*U - cons2*A2*U.*U+cons1*G+cons2*G2;
U_new = reshape(U_new, [m 1]);


end

