function [U_new, k] = adapted_mesh_solver(U, gL, gR, m, eps, h)
%FTBS_SOLVER our ftbs solver
%Solves problems of the form u_t(x,t) + a*u_x(x, t) = 0
%   * U: [m-1, 1] previous iteration
%   * gL: left boundary condition
%   * gR: right boundary condition
%   * m: size of vector
%   * esp: parameter a (determines which sided variant will be used)
%   * h : [m, 1] vector of the distant
%
%   * U_new: new iteration
%   * k: step taken

k = min(h)^2 / (2*eps);

e = ones(m-1, 1);
z = zeros(m-1, 1);
H = sparse(1:m-1, 1:m-1, (h(1:end-1)+h(2:end))/2);

A1 = spdiags([z -e e], [-1 0 1], m-1, m-1);
A2 = spdiags([-e e z], [-1 0 1], m-1, m-1);
U = reshape(U, [m-1 1]);

%estimation of the first and second derivative of U wtr x
dxU = (H^(-1)*A1*U.*U+H^(-1)*A2*U.*U)/2;
dx2U =  H^(-2)*A1*U-H^(-2)*A2*U;

G = zeros(m-1, 1);
G(1) = -k*gL^2/H(1)+k*eps*gL/H(1)^2;
G(end) = k*gR^2/H(end)+k*eps*gR/H(end)^2;

U_new = U + k*eps*dx2U-k*dxU+G;
U_new = reshape(U_new, [m-1 1]);
end
