%% Exercise 1(a)-(b)
clc; clear;

%first stencil
a = -4;
b = 0;
q = 2; %order of derivative
stencil = a:1:b;

%coefficients
coeff_1 = finitedifferences(stencil,q);
%order of accuracy 
C_1 = orderAccuracy(coeff_1,a,b,q);

%second stencil
a = -2;
b = 2;
q = 2;
stencil = a:1:b;

%coefficients
coeff_2 = finitedifferences(stencil,q);
%order of accuracy
C_2 = orderAccuracy(coeff_2,a,b,q);

%% Exercise 1(c)-(d)
clc; clear; close all;
% Calculate optimal coefficients using fdcoeffF and perform convergence

x = -2:1:2; %vector of grid points
u = exp(cos(x)); %vector of the actual function evaluated in the grid points

k = 2; %order of approximated derivative
xbar = 0; 

c = fdcoeffF(k,xbar,x); %coefficents  of the approximation
Order = orderAccuracy(c,-2,2,2);
u_2 = -(cos(xbar).*exp(cos(xbar))+(sin(xbar)).^2.*exp(cos(xbar))); %second derivative of the actual function
u_app = c*u'; %approximation of the second derivative

%plots 
figure
plot(xbar,u_2,'o','Color','b');
hold on;
plot(xbar,u_app,'o','Color','r');
grid on;
legend('u_2', 'u_a_p_p')
title('Exact solution - Approximation')

err = abs(u_2-u_app);
%%
h = zeros(1,7); %vector of different step sizes
for n = 2:8
    h(n-1) = 1/2^n;
end 

%construction of the vector containing the approximation with different
%values of steps
results = zeros(length(h),1);  
for i = 1:length(h)
    x = (xbar-2*h(i)):h(i):(xbar+2*h(i)); %stencil
    c = fdcoeffF(k,xbar,x);
    u = exp(cos(x));
    results(i) = c*u';
end 

figure
loglog(h(:), results(:));
hold on;
loglog(h(:), results(:), 'o');

n_iterate = 1:length(h);
offset_x = 1e-4; 
offset_y = 1e-4; 
for i = 1:length(h)
    text(h(i)+offset_x, results(i)+offset_y, ['n=', num2str(n_iterate(i))]); 
end

title('Approximation as a function of h')
xlabel('h'); ylabel('u_a_p_p');
grid on;

%total error
error = zeros(length(h),1);
for i = 1:length(h)
    error(i) = results(i)+exp(1);
end 

figure
loglog(h(:), error(:));
hold on;
loglog(h(:), error(:), 'o');
title('Total error as a function of h');
xlabel('h'); ylabel('Total error');
grid on;

%%
clc; clear;

i = 1/2;

%first stencil
stencil = -1/2:1:1/2;
coeff = finitedifferences(stencil,0);
C = orderAccuracy(coeff,-i,i,0);
%%
%second stencil
stencil = -3/2:1:3/2;
coeff = finitedifferences(stencil,0);
C = orderAccuracy(coeff,-i,i,0);

%%
%forward stencil
stencil = -1/2:1:3/2;
coeff = finitedifferences(stencil,0);
C = orderAccuracy(coeff,-1/2,3/2,0);

%%
%backward stencil
stencil = -3/2:1:1/2;
coeff = finitedifferences(stencil,0);
C = orderAccuracy(coeff,-3/2,1/2,0);
