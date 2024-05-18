function T_star = picard_lindeloef(t0, t1, y0, fun)
% PICARD_LINDELOEF apply Picard Lindeloef test using fairly arbitrary
% parameters. This function will return an interval where a solution is
% believed to exist, however it does not check that the input function has
% sufficient smoothness to validate that.
%    * T0 beginning of the interval
%    * T1 maximum T value
%    * Y0 initial condition
%    * FUN function
%
%    * T_STAR output T value for which a solution is believed to exist

T_star = 0;
a = 0.1;
while true
% pick a

% optimize S
S = max(abs(fun(0, linspace(y0-a, y0+a, 100))));

T_star_new = a/S;
if (T_star_new < T_star) || (t0 + T_star_new >= t1)
    break
else
    T_star = T_star_new;
    a = 1 + 0.1;
end
end

T_star = min(t1, t0 + T_star);

end
