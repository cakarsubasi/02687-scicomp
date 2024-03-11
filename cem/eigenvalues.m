%% Determine best omega for relaxed Jacobi smoothing
clear;
clc;

m = 5000;
h = 1 / (m+1);
p = m/2:m;
q = (m/2:m)';
lambda_p_q = @(p, q, h) ((cos(h*pi*p) - 1) + (cos(h*pi*q) - 1));

lambda_p_q_2 = @(p, q, h, omega) 1 + omega * 0.5 * lambda_p_q(p, q, h);

lam = lambda_p_q_2(p, q, h, 0.66);
omegas_n = 50;
omegas = linspace(0, 2, omegas_n)';

surf(p, q, lam);

maxes = zeros(omegas_n, 1);
for i = 1:omegas_n
    omega = omegas(i);
    maxes(i) = max(abs(lambda_p_q_2(p, q, h, omega)), [], "all");
end

plot(omegas, maxes);
disp(max(lam, [], "all"));
%%
%lam2 = lambda_p_q_2(omegas);
%surf(lam2);
surf(p, q, lam);
xlabel("p");
ylabel("q");
