function rhs = form_rhs(m,f,u)
%FORM_RHS Summary of this function goes here
%   m: number of points per axis
%   f: function handle
%   u: boundary conditions (m x m matrix)
[rhs, ~, ~] = makerhs(f, "5-point", [0 1 0 1], m, u);
end

