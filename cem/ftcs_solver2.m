function [U_new, k] = ftcs_solver2(U, G, m, eps)
%FTCS_SOLVER our ftcs solver
%   * U: u values for previous iteration [1 x m]
%   * G: boundary values [2 x 1]
%   * m: number of values
%   * eps:

h = 2/(m+1);
k = h^2 / (2*eps);
% pick h and k values
cons = eps*k/h^2;

U = reshape(U, [m, 1]);
G = [G(1); zeros([m, 1]); G(2)];

e = ones(m, 1);

A = spdiags([1*e -2*e 1*e], [-1 0 1], m, m)*cons;
A = [zeros([1, m]); A; zeros([1, m])];

U_new = [0; U; 0] + A*U + G*cons;
U_new = U_new(2:end-1);

end