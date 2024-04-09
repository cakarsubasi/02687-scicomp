% stencil: (a ... b), n: derivative order (integer)
function [coeff] = finitedifferences(stencil,n)
% Compute coefficients for finite difference approximation for the
% derivative of order n based on grid values at points in stencil.
%
% INPUT
% "stencil"  vector [1xk] containing the grid points   
% "n"        order of the derivative that needs to be approximated
% 
% OUTPUT
% "coeff"    vector [kx1] containing the coefficients of the approximated function
%

k = length(stencil);

A = zeros(k);
fact = 1;

for i = 1:k
    for j = 1:k
        A(i,j) = (stencil(j))^(i-1)/fact;
    end 
    fact = fact*i;
end
f = zeros(k,1);
f(n+1) = 1;
coeff = A\f;

end 