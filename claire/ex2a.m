a=0;
b=1;
alpha = -1;
beta = 1.5;
epsilon = 0.1;
tol = 0.0001;
kmax = 50;
n = 30;
h = (b-a)/(n-1);
omega0 = (a-b-alpha+beta)/2;
x = (a:h:b)';
x
xbar =(a+b-alpha-beta)/2*ones(n, 1);
xbar
u0 = x - xbar + omega0*tanh(omega0*(x-xbar)/(2*epsilon));
u0
u = NewtonMethod(u0, epsilon, h, alpha, beta, tol, kmax);
u
plot(x, u)
