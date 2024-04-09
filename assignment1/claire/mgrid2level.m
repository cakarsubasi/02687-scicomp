% exact solution
u_exact = @(x, y) cos(4*pi*x.*y) + sin(4*pi*(x + y));
% u_xx + u_yy (right hand side function)
f = @(x, y) - 32*pi^2*sin(4*pi*(x + y)) - 16*x.^2*pi^2.*cos(4*pi*x.*y) - 16*y.^2*pi^2.*cos(4*pi*x.*y);

m = 25;
h=1/(m+1);
omega = 2/3;


bc = [0 1 0 1];

range = h:h:(1-h);
[Xint,Yint]=meshgrid(range, range);

[F, u_solution, h] = makerhs(f, "5-point", bc, m, u_exact);
Uhat = u_solution(2:(m+1), 2:(m+1));

U = zeros(m^2, 1);

for i=1:10,
  U = smooth(U, omega, m, F);
  U2 = reshape(U, [m, m]);
  E2 = U2-Uhat;
  subplot(1,2,1),
  surf(Xint,Yint,U2);
  set(gca,'fontsize',16);
  xlabel('x');
  ylabel('U');
  title(sprintf('Iter=%4d', i));
  subplot(1,2,2),
  surf(Xint,Yint,E2);
  set(gca,'fontsize',16);
  xlabel('x');
  ylabel('E');
  title(sprintf('Iter=%4d', i));
  set(gcf,'color',[1,1,1]);
  pause(1);
end
pause(5);
% calculate residual
r = F + Amult(U, m);
% coarsen
m_coarse = (m-1)/2;
h_coarse = 1/(m_coarse+1);
r_coarse = coarsen(r, m);
% solve the coarse problem
A_coarse= poisson5(m_coarse);
e_coarse = A_coarse\r_coarse;
% project back on the fine grid
e = interpolate(e_coarse, m);
U = U+e;
U2 = reshape(U, [m, m]);
E2=U2-Uhat;
subplot(1,2,1),
surf(Xint,Yint,U2);
set(gca,'fontsize',16);
xlabel('x');
ylabel('U');
title('After coarse grid projection');
subplot(1,2,2),
surf(Xint,Yint,E2);
set(gca,'fontsize',16);
xlabel('x');
ylabel('E');
title('After coarse grid projection');
set(gcf,'color',[1,1,1]);
pause(5);
% smooth the error again
for i=1:10,
  U = smooth(U, omega, m, F);
  U2 = reshape(U, [m, m]);
  E2 = U2-Uhat;
  subplot(1,2,1),
  surf(Xint,Yint,U2);
  set(gca,'fontsize',16);
  xlabel('x');
  ylabel('U');
  title(sprintf('Iter=%4d', i));
  subplot(1,2,2),
  surf(Xint,Yint,E2);
  set(gca,'fontsize',16);
  xlabel('x');
  ylabel('E');
  title(sprintf('Iter=%4d', i));
  set(gcf,'color',[1,1,1]);
  pause(1);
end
