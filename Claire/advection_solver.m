function [U_new, k] = advection_solver(U, gL, gR, m, eps)
%FTBS_SOLVER our ftbs solver
%Solves problems of the form u_t(x,t) + a*u_x(x, t) = 0
%   * U: [m-1, 1] previous iteration
%   * gL: left boundary condition
%   * gR: right boundary condition
%   * m: size of vector
%   * a: parameter a (determines which sided variant will be used)
%
%   * U_new: new iteration
%   * k: step taken

h = 2/(m+1);
k = h^2 / (2*eps);

e = ones(m-1, 1);
z = zeros(m-1, 1);
A1 = spdiags([1*e -2*e 1*e], [-1 0 1], m-1, m-1);
A2 = spdiags([-1*e z 1*e], [-1 0 1], m-1, m-1);
cons1 = k*eps/h^2;
cons2 = k/(2*h);
U = reshape(U, [m-1 1]);
G = zeros(m-1, 1);
G(1) = cons1*gL+cons2*gL*gL;
G(end) = cons1*gR-cons2*gR*gR;
U_new = U + cons1*A1*U - cons2*A2*U.*U+G;
U_new = reshape(U_new, [m-1 1]);



end
