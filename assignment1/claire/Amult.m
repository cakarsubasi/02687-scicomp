function AU = Amult(U, m)
  ##%AMULT Summary of this function goes here
  ##%   For exercise 3.1
  ##handle = poisson5_f(m);
  ##AU = -pcg(handle, U, 1e-6, 200);
  a = @(x) -4*x + [0; x(1:end-1)] + [x(2:end); 0];
  b = @(x) x;

  ax_i = @(x, i) a(x(i*m+1:(i+1)*m));
  bx_i = @(x, i) b(x(i*m+1:(i+1)*m));

  top = @(x) ax_i(x, 0) + bx_i(x, 1); % n^2 -> n
  bot = @(x) ax_i(x, m-1) + bx_i(x, m-2); % n^2 -> n
  centers = @(x, i) ax_i(x, i) + bx_i(x, i-1) + bx_i(x, i+1); % n^2 -> n

  AU = zeros(m*m, 1);
  AU(1:m) = top(U);
  AU(m*(m-1)+1:end) = bot(U);
  for i = 1:m-2
    AU(m*i+1:m*(i+1)) = centers(U, i);
  end
  AU = -(m+1)^2*AU;
end
