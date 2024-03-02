function AU = Amult(U, m)
%AMULT Summary of this function goes here
%   For exercise 3.1
handle = poisson5_f(m);

AU = -pcg(handle, U, 1e-6, 200);
end

