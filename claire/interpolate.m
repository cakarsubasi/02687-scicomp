function R=interpolate(Rc,m)
  R = zeros(m^2, 1);
  mc = (m-1)/2;
  R(1) = Rc(1);
  R(m) = Rc(mc);
  R(m^2-m+1)  = Rc(mc^2-mc+1);
  R(m^2) = Rc(mc^2);
  for i=2:(m-1)
    ic = floor(i/2);
    if mod(i, 2) == 0
      R(i) = Rc(ic);
      R(m^2-m+i) = Rc(mc^2-mc+ic);
      R(m*i) = Rc(mc*ic);
      R(m*(i-1)+1) = Rc(mc*(ic-1)+1);
    else
      R(i) = 1/2*(Rc(ic)+Rc(ic+1));
      R(m^2-m+i) = 1/2*(Rc(mc^2-mc+ic)+Rc(mc^2-mc+ic+1));
      R(m*i) = 1/2*(Rc(mc*ic)+Rc(mc*(ic+1)));
      R(m*(i-1)+1) = 1/2*(Rc(mc*ic+1)+Rc(mc*(ic-1)+1));
    endif
  endfor
  for i=2:(m-1)
    for j=2:(m-1)
      ic = floor(i/2);
      jc = floor(j/2);
      k = m*(i-1)+j;
      if mod(i, 2) == 1
        if mod(j, 2) == 1
          R(k) = 1/4*(Rc(mc*(ic-1)+jc)+Rc(mc*(ic-1)+jc+1)+ Rc(mc*ic+jc)+ Rc(mc*ic+jc+1));
         else
          R(k) = 1/2*(Rc(mc*(ic-1)+jc)+Rc(mc*ic+jc));
        endif
      else
        if mod(j, 2) == 1
          R(k) = 1/2*(Rc(mc*(ic-1)+jc)+Rc(mc*(ic-1)+jc+1));
         else
          R(k) = Rc(mc*(ic-1)+jc);
        endif
      endif
    endfor
  endfor
end
