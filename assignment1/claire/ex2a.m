##constante
a=0;
b=1;
alpha = -1;
beta = 1.5;
epsilon = 0.1;
tol = 0.00001;
kmax = 10000;


##calculation of the solution for u=1
n=100;
h = 1/(n-1);
omega0 = (a-b-alpha+beta)/2;
x = (a:h:b)';
xbar =(a+b-alpha-beta)/2*ones(n, 1);
u0 = x - xbar + omega0*tanh(omega0*(x-xbar)/(2*epsilon));
[u, error, delta] = NewtonMethod(u0, epsilon, h, alpha, beta, tol, kmax);
plot(x, u) %plot the solution
n_delta = length(delta);
semilogy(1:n_delta, delta); %plot the convergence


##calculation of the errors with regard to h
nb_iter= 20;
N=logspace(1, 3.5, nb_iter);
H = zeros(1, nb_iter);
e = ones(1, nb_iter);
for k = 1:nb_iter
  n = round(N(k));
  h = (b-a)/(n-1);
  H(k) = h;
  omega0 = (a-b-alpha+beta)/2;
  x = (a:h:b)';
  xbar =(a+b-alpha-beta)/2*ones(n, 1);
  u0 = x - xbar + omega0*tanh(omega0*(x-xbar)/(2*epsilon));
  [u, err, delta] = NewtonMethod(u0, epsilon, h, alpha, beta, tol, kmax);
  e(k) = err;
endfor

p = polyfit(H, e, 1);
logH = log(H);
[p, error] = polyfit(logH(4:nb_iter), log(e(4:nb_iter)), 1);
P = polyval(p, logH);
loglog(H, e)
hold on
loglog(H, exp(P))
xlabel('log(h)')
ylabel('log(E)')
legend('Global error estimate', 'Linear regression','Location','southeast')
