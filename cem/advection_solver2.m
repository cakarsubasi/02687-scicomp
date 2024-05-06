function [U_new, k] = advection_solver2(U, G, m, eps, t)
%ADVECTION_SOLVER2 
% Solves the advection problem in 4.3
%   * U: [m, 1] previous iteration
%   * G: [m, 1] boundary conditions, only the first or last element is used
%   * m: size of vector
%   * eps: parameter eps
%   * t: current time
%
%   * U_new: new iteration
%   * k: step taken

h = 2/(m+1);
k = h^2 / (2*eps);
if (t + k) >= 1.6037/pi
    k = 1.6037/pi - t;
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

