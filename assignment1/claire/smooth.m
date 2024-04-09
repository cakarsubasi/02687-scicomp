function Unew=smooth(U,omega,m,F)
  U = reshape(U, [m, m]);
  Unew = [zeros(1,m+2) ; [zeros(m,1), U, zeros(m,1)] ; zeros(1,m+2)];
  Fnew = reshape(F, m, m);
  I = 2:(m+1);
  J = 2:(m+1);
  h = 1/(m+1);
  Unew(I,J) = (1-omega)*Unew(I, J) + 0.25*omega*(Unew(I-1,J) + Unew(I+1,J) + Unew(I,J-1) + Unew(I,J+1) - h^2 * Fnew(I-1,J-1));
  Unew = reshape(Unew(I, J), m^2, 1);
end
