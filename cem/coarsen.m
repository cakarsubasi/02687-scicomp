function Rcoarse = coarsen(R, m)
%COARSEN Summary of this function goes here
%   Detailed explanation goes here

% coarsen
mc = (m-1)/2;

Rcoarse = reshape(R, m, m);
Rcoarse = Rcoarse(2:2:end, 2:2:end);
assert(all(size(Rcoarse)==[mc, mc]));
Rcoarse = reshape(Rcoarse, mc*mc, 1);

end

