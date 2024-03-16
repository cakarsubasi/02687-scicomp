function [u, err, l_delta] = NewtonMethod(u0, eps, h, alpha, beta, tol, k_max)

  #Initialization of u
  u = u0;
  n = size(u0)(1);
  G = zeros(n, 1);
  large_u = [0; u; 0];
  l_delta=zeros(100, 1);
  k=0;

  #Calculation of G
  for j=2:(n-1)
    G(j) = eps*(u(j-1)-2*u(j)+u(j+1))/(h^2)+u(j)*((u(j+1)-u(j-1))/(2*h)-1);
  endfor


  while norm(G)>tol && k<k_max

    #Calculation of J
    d1 = eps*ones(n-1, 1)/(h^2)-u(2:n)/(2*h);
    d2 = (-2*eps/(h^2)-1)*ones(n, 1)+(large_u(3:(n+2))-large_u(1:n))/(2*h);
    d3 = eps*ones(n-1, 1)/(h^2)+u(1:(n-1))/(2*h);
    d4 = [0; d3];
    J = spdiags(d1, -1, n, n)+spdiags(d2, 0, n, n)+spdiags(d4, 1, n, n);

    #Calculation of delta
    delta = J\G;

    #Correction of u
    u = u-delta;
    u(1)=alpha;
    u(n)=beta;

    #update the variable
    k = k+1;
    for j=2:(n-1)
      G(j) = eps*(u(j-1)-2*u(j)+u(j+1))/(h^2)+u(j)*((u(j+1)-u(j-1))/(2*h)-1);
    endfor
    large_u = [0; u; 0];
    l_delta(k) = norm(delta, inf);
  endwhile
  err = norm(J\G, Inf);
