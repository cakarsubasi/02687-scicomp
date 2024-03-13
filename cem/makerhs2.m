function [rhs, u_sol, h] = makerhs2(func, bc, m, bc_matrix_or_solution, del_f)
% for the 9-point case where del of f is analytically known
a = bc(1); b = bc(2); c = bc(3); d = bc(4); 

h = (b-a)/(m+1);
% warning: this only supports uniform bcs

X = linspace(a, b, m+2); % all points
Y = linspace(c, d, m+2);
[X,Y] = meshgrid(X,Y);
Xindices = 2:m+1; % interior points
Yindices = 2:m+1;
Xint = X(Xindices,Yindices);       % interior points
Yint = Y(Xindices,Yindices);

rhs = func(Xint, Yint);

% convert from function handle if needed
if isa(bc_matrix_or_solution, 'function_handle')
    bc_matrix = bc_matrix_or_solution(X, Y);
else
    bc_matrix = bc_matrix_or_solution;
end
u_sol = bc_matrix;

rhs_err = (h^2 / 12) * del_f(Xint, Yint);
rhs = rhs + rhs_err;

% adjust the rhs to include boundary terms:
rhs(:,1) = rhs(:,1) - (4*bc_matrix(Xindices, 1) + ...
    bc_matrix(Xindices-1, 1) + bc_matrix(Xindices+1, 1))/6/h^2;  % top
rhs(:,m) = rhs(:,m) - (4*bc_matrix(Xindices,m+2) + ...
    bc_matrix(Xindices-1,m+2) + bc_matrix(Xindices+1,m+2))/6/h^2; % bottom
rhs(1,:) = rhs(1,:) - (4*bc_matrix(1,Yindices) + ...
    bc_matrix(1,Yindices-1) + bc_matrix(1,Yindices+1))/6/h^2;   % left
rhs(m,:) = rhs(m,:) - (4*bc_matrix(m+2,Yindices) + ...
    bc_matrix(m+2,Yindices-1) + bc_matrix(m+2,Yindices+1))/6/h^2; % right

% fix the corners
rhs(1,1)      = rhs(1, 1)     + bc_matrix(1, 1)    /6/h^2; % top-left
rhs(1,end)    = rhs(1, end)   + bc_matrix(1, end)  /6/h^2; % bottom-left
rhs(end,1)    = rhs(end, 1)   + bc_matrix(end, 1)  /6/h^2; % top-right
rhs(end, end) = rhs(end, end) + bc_matrix(end, end)/6/h^2; % bottom-right

% convert matrix into column vector
rhs = reshape(rhs, m*m, 1);
end
