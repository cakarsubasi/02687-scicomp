function Einterp = interpolate(Ecoarse, mc)
%INTERPOLATE Summary of this function goes here
%   Detailed explanation goes here

m = mc*2 + 1;
Einterp=zeros(m, m);
Einterp(2:2:end, 2:2:end) = reshape(Ecoarse, mc, mc);

%%interp2(Xq, Yq, reshape(coarsen(U_sol, m), mc, mc), X, Y);

for i=1:2:m
  if(i>1)
    e_left = Einterp(:, i-1);
  else
    e_left = Einterp(:, 2);
  end
  if(i<m)
    e_right = Einterp(:, i+1);
  else
    e_left = Einterp(:, m-1);
  end
  Einterp(:, i) = (e_left+e_right)/2;
end

for i=1:2:m
  if(i>1)
    e_top = Einterp(i-1, :);
  else
    e_top = Einterp(2, :);
  end
  if(i<m)
    e_bottom = Einterp(i+1, :);
  else
    e_bottom = Einterp(m-1, :);
  end
  Einterp(i, :) = (e_top+e_bottom)/2;
end
Einterp = reshape(Einterp, m*m, 1);

end

