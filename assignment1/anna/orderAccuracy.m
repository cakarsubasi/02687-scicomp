function Order = orderAccuracy(coeff,a,b,q)
% Compute the order of accuracy of a stencil  
%
% INPUT
% "coeff"  vector [kx1] containing the approximating coefficents   
% "q"      order of the approximated derivative 
% "a"      left endpoint of the stencil
% "b"      right endpoint of the stencil
% 
% OUTPUT
% "Order"  order of accuracy of the stencil
%

C = 0;

k = 0;
for i = a:1:b
    k = k+1;
end 

n = q+1; %starting the computation from q+1 as for n<q C is certainly null
tol = 1e-3; %threshold to better estimate when C is zero

while C < tol && C > -tol
    j = 0;
    for i = a:1:b
        j = j+1;
        C = C + coeff(j)*i^n;
    end 
    n = n+1;
end 

%order of accuracy 
Order = (n-1)-q;

end 