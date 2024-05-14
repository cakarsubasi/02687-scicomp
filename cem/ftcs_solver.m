function [U_new, k] = ftcs_solver(U, G, m, eps, step_fac, kmax)
%FTCS_SOLVER our ftcs solver
%   * U: u values for previous iteration [1 x m]
%   * G: boundary values [2 x 1]
%   * m: number of values
%   * eps:
%   * step_fac: limit or increase allowed steps, setting it above 1.0 will
%   cause instability
%   * kmax: maximum allowed step, use it to stop an iteration at a specific
%   point

% pick h and k values
h = 2/(m+1);
k = h^2 / (2*eps)*step_fac;
if k > kmax
    k = kmax;
end

cons = eps*k/h^2;

U = reshape(U, [m, 1]);
G = [G(1); zeros([m-2, 1]); G(end)];


e = ones(m, 1);

A = spdiags([1*e -2*e 1*e], [-1 0 1], m, m);

U_new = U + A*U*cons + G*cons;
end

