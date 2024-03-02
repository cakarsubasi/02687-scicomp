function func = poisson5_f(m)
%POISSON5_F 5-point poisson represented as a computational graph
%   Detailed explanation goes here
% We are trying to compute a matrix multiplication that looks like this:
% [ A B 0 0 0 ] [x1]
% [ B A B 0 0 ] [...]
% [ 0 B A B 0 ] [...]
% [ 0 0 B A B ] [...]
% [ 0 0 0 B A ] [x_m^2]
% (m^2 x m^2) * (m^2 x 1)

a = @(x) -4*x + [0; x(1:end-1)] + [x(2:end); 0];
b = @(x) x;

ax_i = @(x, i) a(x(i*m+1:(i+1)*m));
bx_i = @(x, i) b(x(i*m+1:(i+1)*m));

top = @(x) ax_i(x, 0) + bx_i(x, 1); % n^2 -> n
bot = @(x) ax_i(x, m-1) + bx_i(x, m-2); % n^2 -> n
centers = @(x, i) ax_i(x, i) + bx_i(x, i-1) + bx_i(x, i+1); % n^2 -> n

    function c = calculate(x)
        c = zeros(m*m, 1);
        c(1:m) = top(x);
        c(m*(m-1)+1:end) = bot(x);
        for i = 1:m-2
            c(m*i+1:m*(i+1)) = centers(x, i);
        end
    end

func = @(x) -((m+1)^2) * calculate(x);
end