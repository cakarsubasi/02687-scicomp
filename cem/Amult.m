function AU = Amult(U, m)
%AMULT Left multiply U with A without a matrix
%   For exercise 3.1
handle = poisson5_f(m);
AU = handle(U);
end

