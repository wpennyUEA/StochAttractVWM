function [h] = mda_diag_hess (x,mda,load)
% Numerical approx to diagonal Hessian using central differences
% FORMAT [h] = mda_diag_hess (x,mda,load)

P = length(x);
dx = 0.001;
for p=1:P,
    [x1,x2] = exp_points(x,p,dx);
    g1 = grad(x1,mda,load,p,dx);
    g2 = grad(x2,mda,load,p,dx);
    h(p) = (g1-g2)/(2*dx);
end

end    

%---------------------------------------------
function [g] = grad (x,mda,load,p,dx)

[x1,x2] = exp_points(x,p,dx);
E1 = mda_neglogjoint(x1,mda,load);
E2 = mda_neglogjoint(x2,mda,load);
g = (E1-E2)/(2*dx);

end

%---------------------------------------------
function [x1,x2] = exp_points(x,p,dx)

x1 = x;
x1(p) = x(p)+dx;
x2 = x;
x2(p) = x2(p)-dx;

end