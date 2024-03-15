function Unew = jacobi_iter(U, omega, m, F)
%JACOBIS Summary of this function goes here
% Perform Jacobi iteration
%   Detailed explanation goes here
% TODO. For exercise 3.2
% need to fix the boundries

Iint = 2:(m+1);
Jint = 2:(m+1);
h = 1/(m+1);

% convert from function handle if needed
if isa(F, 'function_handle')
    range = 0:h:1; 
    F = F(range, range');
else
end

Unew = zeros(m+2, m+2);
Uold = zeros(m+2, m+2);
Uold(2:m+1, 2:m+1) = reshape(U, m, m);

%Unew(:,1) = F(:, 1);
%Unew(:,m+2) = F(:, m+2);
%Unew(1,:) = F(1,:);
%Unew(m+2,:) = F(m+2,:);
Unew(Iint, Jint) = (1-omega)*Uold(Iint, Jint) + ...
    omega* (0.25*(Uold(Iint-1,Jint) + Uold(Iint+1,Jint) + Uold(Iint,Jint-1) + Uold(Iint,Jint+1) ...
    - h^2 * F(Iint,Jint)));

Unew = reshape(Unew(Iint, Jint), m*m, 1);

end

