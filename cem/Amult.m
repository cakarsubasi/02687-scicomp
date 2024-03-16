function AU = Amult(U, m)
%AMULT Matrix-free Poisson solver
%   For exercise 3.1
handle = poisson5_f(m);
AU = handle(U);
%AU = -pcg(handle, U, 1e-6, 200);
end

