%% Example



%% exercise 1(c)
% symbolic solver:
syms x
f = exp(cos(x));
f_prime_prime = diff(f, 2);
% exp(cos(x))sin(x)^2 - exp(cos(x))*cos(x)
% In this case, we know the exact derivative,
% so we can lazily calculate the truncation error
actual = double(subs(f_prime_prime, x, 0));
h = 1e-3;

xbars = [0 -1 -2 -3 -4];
xbars2 = -1/2:1:3/2; % 3
xbars3 = -3/2:1:5/2; % 4
coeff = fdcoeffF(2, 0, xbars);
coeff2 = fdcoeffF(2, 0, xbars2);
coeff3 = fdcoeffF(2, 0, xbars3);

est = estimator(@fun, xbars, coeff, h);
act2 = exp(cos(0))*sin(0)^2 - exp(cos(0))*cos(0);

h_values = (2*ones(1, 10)).^(-linspace(1, 10, 10));

estimates = arrayfun(@(h) estimator(@fun, xbars, coeff, h), ...
    h_values);
estimates2 = arrayfun(@(h) estimator(@fun, xbars2, coeff2, h), ...
    h_values);
estimates3 = arrayfun(@(h) estimator(@fun, xbars3, coeff3, h), ...
    h_values);

disp(estimates2)

%% plots
loglog(h_values, (act2 - estimates).^2, "o-")
hold on;
loglog(h_values, (act2 - estimates2).^2, "o-")
loglog(h_values, (act2 - estimates3).^2, "o-")
xlabel("h")
ylabel("Total error (\tau + \epsilon)")

%%
function c = fun(x)
c = exp(cos(x));
end

function estimate = estimator(func, xbar, coeff, h)
    estimate = sum(coeff .* func(xbar*h))/h^2;
end

