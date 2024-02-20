function n = C(alpha, beta, q, A)

  n=q;
  Cn = 0;
  while norm(Cn)<0.0001 && n<10
    if n==q
      Cn = -1;
    else
      Cn=0;
    endif
    facto = factorial(n);
    for m=alpha:beta
      Cn = Cn + A(m-alpha+1)*(m^n)/facto;
    endfor
    n=n+1;
  end
  n=n-1;

