function Rcoarse = coarsen(R, m)
%COARSEN Coarsen a given RHS vector by one order of magnitude

% coarsen
mc = (m-1)/2;

Rcoarse = reshape(R, m, m);
Rcoarse = Rcoarse(2:2:end, 2:2:end);
assert(all(size(Rcoarse)==[mc, mc]));
Rcoarse = reshape(Rcoarse, mc*mc, 1);

end

