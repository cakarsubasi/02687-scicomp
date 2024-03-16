function Rc=coarsen2(R,m)
% Claire's Coarsen
  mc = (m-1)/2;
  hc = 1/(mc-1);
  Rc = zeros(mc^2, 1);
  for ic = 1:mc
    for jc = 1:mc
      kc =(ic-1)*mc+jc;
      k = (2*ic-1)*m+2*jc;
      Rc(kc) = R(k-1)/8+R(k+1)/8+R(k)/4+R(k+m)/8+R(k-m)/8+R(k-m-1)/16+R(k-m+1)/16+R(k+m-1)/16+R(k+m+1)/16;
    end
  end
end

