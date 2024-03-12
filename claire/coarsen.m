function Rc=coarsen(R,m)
  mc = (m+1)/2;
  hc = 1/(mc+1);
  Rc = zeros(mc^2, 1);
  Rc(1) = R(1)/4+R(2)/8+R(m+1)/8+R(m+2)/16;
  Rc(mc) = R(m)/4+R(m-1)/8+R(2*m)/8+R(2*m-1)/16;
  Rc(mc^2-mc+1) = R(m^2-m+1)/4+R(m^2-m+2)/8+R(m^2-2*m+1)/8+R(m^2-2*m+2)/16;
  Rc(mc^2) = R(m^2)/4+R(m^2-1)/8+R(m^2-m)/8+R(m^2-m-1)/16;
  for ic = 2:(mc-1)
    i = 2*ic-1;
    Rc(ic) = R(i-1)/8+R(i)/4+R(i+1)/8+R(i+m)/8+R(i+m-1)/16+R(i+m+1)/16;
    Rc(mc^2-mc+ic) = R(m^2-m+i-1)/8+R(m^2-m+i)/4+R(m^2-m+i+1)/8+R(m^2-2*m+i)/8+R(m^2-2*m+i-1)/16+R(m^2-2*m+i+1)/16;
    Rc(mc*(ic-1)+1) = R(m*(i-2)+1)/8+R(m*(i-1)+1)/4+R(m*i+1)/8+R(m*(i-1)+2)/8+R(m*(i-2)+2)/16+R(m*i+2)/16;
    Rc(mc*ic) = R(m*(i-1))/8+R(m*i)/4+R(m*(i+1))/8+R(m*i-1)/8+R(m*(i-1)-1)/16+R(m*(i+1)-1)/16;
    for jc = 2:(mc-1)
      kc =(ic-1)*mc+jc;
      k = (i-1)*m+2*jc-1;
      Rc(kc) = R(k-1)/8+R(k+1)/8+R(k)/4+R(k+m)/8+R(k-m)/8+R(k-m-1)/16+R(k-m+1)/16+R(k+m-1)/16+R(k+m+1)/16;
    endfor
  endfor
end

