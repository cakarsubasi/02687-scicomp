function Unew = jacobis(U, omega, m, F)
%JACOBIS Summary of this function goes here
%   Detailed explanation goes here
% TODO. For exercise 3.2
% need to fix the boundries
[r, c] = size(U);
Unew = [zeros(1,c+2) ; [zeros(r,1), U, zeros(r,1)] ; zeros(1,c+2)];

maxiter = 200;
I = 2:(m+1);
J = 2:(m+1);
h = 1/(m+1);
for iter=0:maxiter
Unew(I,J) = (1-omega)*Unew(I, J) + omega*(0.25*(Unew(I-1,J) + Unew(I+1,J) + Unew(I,J-1) ...
+ Unew(I,J+1) - h^2 * F(I,J)));
end
Unew = Unew(I,J);
end

