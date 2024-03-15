%% exact solution and RHS
%u=@(x,y) exp(pi*x).*sin(pi*y)+0.5*(x.*y).^2;
%f=@(x,y) x.^2+y.^2;
u = @(x, y) cos(4*pi*x.*y) + sin(4*pi*(x + y));
f = @(x, y) - 32*pi^2*sin(4*pi*(x + y)) - 16*x.^2*pi^2.*cos(4*pi*x.*y) ...
    - 16*y.^2*pi^2.*cos(4*pi*x.*y);
m=2^6-1;
U =zeros(m*m,1);
F = form_rhs(m,f,u); 
epsilon = 1.0E-10;
omega = 0.66;

h = 1/(m-1);
[X,Y]=meshgrid(0:h:1);
U_sol = u(X, Y);
mc=(m-1)/2;
hc = 1/(mc-1);
U_sol_r2 = reshape(interpolate(coarsen(U_sol, m), mc), m, m);
plotU(m, U_sol_r2);
%%

for i=1:100
    R =F+Amult(U,m);
    fprintf('*** Outer iteration: %3d, rel. resid.: %e\n', ...
        i, norm(R,2)/norm(F,2));
    if(norm(R,2)/norm(F,2) < epsilon)
        break;
    end
    U=Vcycle(U,omega,3,m,F);
    plotU(m,U);
    pause(.5);
end

function Unew=Vcycle(U,omega,nsmooth,m,F)
% Approximately solve: A*U = F
h=1.0/(m+1);
l2m=log2(m+1);
assert(l2m==round(l2m));
assert(length(U)==m*m);
if(m==1)
    % if we are at the coarsest level
    % TODO: solve the only remaining equation directly!
    Unew = 1/F;
else
    % 1. pre-smooth the error
    %    perform <nsmooth> Jacobi iterations
    for i = 1:nsmooth
        U = smooth(U,omega,m,F);
    end
    % 2. calculate the residual
    R = F + Amult(U, m);
    % 3. coarsen the residual
    Rcoarse = coarsen(R, m);
    % 4. recurse to Vcycle on a coarser grid
    mc=(m-1)/2;
    Ecoarse=Vcycle(zeros(mc*mc,1),omega,nsmooth,mc,-Rcoarse);
    % 5. interpolate the error
    Einterp = interpolate(Ecoarse, mc);
    % 6. update the solution given the interpolated error
    % TODO: check the sign
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