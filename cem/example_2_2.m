% solve u''(t) = -sin(u(t)) for 0 < t < 2pi
% u(0) = 0.7, u(2pi) = 0.7
%% preamble
% new comment
% BCs
a = 0;
b = 2*pi;
alpha = 0.7;
beta = 0.7;

m = 50;
h = (b-a)/(m+1);

X = a:h:b;
Xint = a+h:h:b-h; % interior points
e = ones(m, 1);

%% setup
% sparse grid diagonal

theta = zeros(1, m);

% 1/h^2 [-2  1 ...
%       [ 1 -2  1 ...
%       [ 0  1 -2  1 ...
%       [ ...
%% iterate

maxiter = 500;
nplot = 10;
error = zeros(1, maxiter);
for iter=1:maxiter
J = jvector(theta, alpha, beta, h);
%theta(1) = alpha;
%theta(end) = beta;

G = gvector(theta, [alpha beta], h);

update = J\(-G);
update = update';
theta = theta + update;

theta(1) = alpha;
theta(end) = beta;
      
   %error(iter+1) = max(abs(u-ustar));
   if mod(iter,nplot)==0
      % plot u every nplot iterations
      plot(theta)
      hold on;
      title(sprintf("Iteration %4i", iter));
      %title(sprintf('%s: Iteration %4i, error = %9.3e',...
      %       method,iter,error(iter+1)),'FontSize',15)
      drawnow
      pause(.1)
   end
end
%%
plot(Xint, theta)

%% G vector for Newton Method
function G = gvector(theta, bc, h)
    alpha = bc(1); beta = bc(2);
    lhs = [alpha theta beta];
    % use basic 2nd order approximation
    lhs = lhs(1:end-2) - 2*lhs(2:end-1) + lhs(3:end);
    lhs = 1/h^2 * lhs;
    % function is problem dependent
    G = lhs + sin(theta);
    G = G';
end
% Jacobian
function J = jvector(theta, alpha, beta, h)
m = length(theta);
% function is problem specific
j = cos(theta);
j(1) = cos(alpha);
j(end) = cos(beta);
J = h^2 * diag(j);

e = ones(m, 1);
A = spdiags([e -2*e e], [-1 0 1], m, m);

J = 1/h^2 * (J + A);
end