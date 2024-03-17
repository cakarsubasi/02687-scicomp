% solve u''(x) = f(x) for 0 < x < 1
% u(0) = alpha, u(1) = beta

% f(x) is specified
f = @(x) exp(x);
% BCs
a = 0;
b = 1;
alpha = exp(a);
beta = exp(b);

m = 10;
h = (b-a)/(m+1);

X = a+h:h:b-h; % interior points
e = ones(m, 1);

A = spdiags([e -2*e e], [-1 0 1], m, m)/h^2; % sparse grid diagonal
% 1/h^2 [-2  1 ...
%       [ 1 -2  1 ...
%       [ 0  1 -2  1 ...
%       [ ...
F = f(X)'; % first order approximation
% use 1/h^2 (U-1 - 2U + U+1)



F(1) = F(1) - alpha/h^2;
F(end) = F(end) - beta/h^2;

U = A\F;
U = [alpha; U; beta];

plot(U);