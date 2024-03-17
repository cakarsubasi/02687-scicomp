function [U_sol, iter] = multigrid(U, m, F, tol, maxiter, to_plot)
%MULTIGRID Solve the Poisson problem using the multigrid method
%   U: initial guess (size: [m^2, 1])
%   m: number of grid points per axis (needs to be 2^n-1 for positive n
%      integer)
%   F: right hand side
%   tol: tolerance
%   maxiter: maximum number of iterations before giving up
%   to_plot: (bool) whether to plot each iteration

omega = 2/3;
iter = maxiter + 1;
for i=1:maxiter
    R =F+Amult(U,m);
    res = norm(R,2)/norm(F,2);
    fprintf('*** Outer iteration: %3d, rel. resid.: %e\n', ...
        i, res);
    if(res < tol)
        iter = i;
        break;
    end
    U=Vcycle(U,omega,3,m,F);
    if to_plot
        plotU(m,U);
        pause(.5);
    end
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