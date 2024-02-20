### Exercice 1
### a.

A = zeros(5);
facto = 1;
for i=1:5
  for j=1:5
    A(i, j) = ((1-j)^(i-1))/facto;
  endfor
  facto = facto*i;
endfor

D = zeros(5);
facto = 1;
for i=1:5
  for j=1:5
    D(i, j) = ((j-3)^(i-1))/facto;
  endfor
  facto = facto*i;
endfor

B = zeros(5, 1);
B(3)=1;

E = A\B
F=D\B

#b.
C(0, 4, 2, E)
C(-2, 2, 2, F)
#c.
##xbar = 0;
##x = [-1 -0.5 0 0.5 1]
##e = zeros(1, 8)
##for s = 2:9
##  x = cat(2, x(1:s), [-0.5^s 0 0.5^s], x((2*s-s+2):(2*s+1)));
##  U = exp(cos(x));
##  c = fdcoeffF(2, 0, x);
##  u_second =c*U';
##  e(s-1) = abs(u_second+exp(1));
##endfor
##plot(2:9, log(e))

## d.
xbar = 0;
h = 0.5;
x = [-h/2 h/2];
e = zeros(1, 7);

##for s = 2:8
##  x = cat(2, [-s*h+h/2], x, [s*h-h/2]);
##  U = exp(cos(x));
##  y = Lagrange_interpolation(0, x, U)
##  abs(y-exp(1))
##endfor




