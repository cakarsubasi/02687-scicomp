function U_est = jacobi(U, omega, m, F, maxiter, tol)
%JACOBI Summary of this function goes here
%   Detailed explanation goes here
Unew = U;
for i = 1:maxiter
    Unew = jacobi_iter(U, omega, m, F);
end
U_est = Unew;

end

