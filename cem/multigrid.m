function U_sol = multigrid(U, m, F, tol)
%MULTIGRID Summary of this function goes here
%   Detailed explanation goes here

omega = 2/3;

for i=1:100
    R =F+Amult(U,m);
    fprintf('*** Outer iteration: %3d, rel. resid.: %e\n', ...
        i, norm(R,2)/norm(F,2));
    if(norm(R,2)/norm(F,2) < tol)
        break;
    end
    U=Vcycle(U,omega,3,m,F);
    plotU(m,U);
    pause(.5);
end
U_sol = U;

end

function Unew=Vcycle(U,omega,nsmooth,m,F)
% Approximately solve: A*U = F
h=1.0/(m+1);
l2m=log2(m+1);
assert(l2m==round(l2m));
assert(length(U)==m*m);
if(m==1)
    % if we are at the coarsest level
    % solve the only remaining equation directly!
    A = poisson5(m); % this returns one number
    Unew = A\F;
else
    % 1. pre-smooth the error
    %    perform <nsmooth> Jacobi iterations
    for i = 1:nsmooth
        U = smooth(U,omega,m,F);
    end
    % 2. calculate the residual
    R = F + Amult(U, m);
    % 3. coarsen the residual
    Rcoarse = coarsen2(R, m); % Claire's coarsen
    % 4. recurse to Vcycle on a coarser grid
    mc=(m-1)/2;
    Ecoarse=Vcycle(zeros(mc*mc,1),omega,nsmooth,mc,-Rcoarse);
    % 5. interpolate the error
    Einterp = interpolate(Ecoarse, mc);
    % 6. update the solution given the interpolated error
    U = U - Einterp;
    % 7. post-smooth the error
    %    perform <nsmooth> Jacobi iterations
    for i = 1:nsmooth
        U = smooth(U,omega,m,F);
    end
    Unew = U;
end
end

function plotU(m,U)
h=1/(m+1);
x=linspace(h,1-h,m);
y=linspace(h,1-h,m);
[X,Y]=meshgrid(x,y);
surf(X, Y, reshape(U,[m,m])');
shading interp;
title('Computed solution');
xlabel('x');
ylabel('y');
zlabel('U');
end