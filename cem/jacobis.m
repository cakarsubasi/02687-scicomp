function Unew = jacobis(Unew, omega, m, F)
%JACOBIS Summary of this function goes here
% Perform Jacobi iteration
%   Detailed explanation goes here
% TODO. For exercise 3.2
% need to fix the boundries

maxiter = 200;


Iint = 2:(m-1);
Jint = 2:(m-1);
h = 1/(m+1);

for iter=0:maxiter
Unew(Iint, Jint) = (1-omega)*Unew(Iint, Jint) + ...
    omega* (0.25*(Unew(Iint-1,Jint) + Unew(Iint+1,Jint) + Unew(Iint,Jint-1) + Unew(Iint,Jint+1) ...
    - h^2 * F(Iint,Jint)));
end
Unew(:,1) = F(Xindices, 1)/h^2;
Unew(:,m) = F(Xindices, m)/h^2;
Unew(1,:) = F(1,Yindices)/h^2;
Unew(m,:) = F(m+2,Yindices)/h^2;

end

