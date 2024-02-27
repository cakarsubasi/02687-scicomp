
function u = NewtonMethod(u0, eps, h, alpha, beta, tol, k_max)
  u = u0;
  n = size(u0, 1);
  G = zeros(n, 1);
  u_prime = [0; u; 0];
  for j=2:(n-1)
    G(j) = eps*(u(j-1)-2*u(j)+u(j+1))/(h^2)+u(j)*((u(j+1)-u(j-1))/(2*h)-1);
  end
  k=0;
  G
  norm(G)
  while norm(G)>tol && k<k_max
    d1 = eps*ones(n-1, 1)/(h^2)-u(2:n)/(2*h);
    d2 = (-2*eps/(h^2)-1)*ones(n, 1)+(u_prime(3:(n+2))-u_prime(1:n))/(2*h);
    d3 = eps*ones(n-1, 1)/(h^2)+u(1:(n-1))/(2*h);
    J = diag(d1, -1)+diag(d2, 0)+diag(d3, 1);
    delta = J\G;
    u = u-delta;
    u(1)=alpha;
    u(n)=beta;
    k = k+1;
    for j=2:(n-1)
      G(j) = eps*(u(j-1)-2*u(j)+u(j+1))/(h^2)+u(j)*((u(j+1)-u(j-1))/(2*h)-1);
    end
    u_prime = [0; u; 0];
  end
  k

