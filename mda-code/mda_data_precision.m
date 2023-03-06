function [h] = mda_data_precision (x,mda,load)
% Compute pointwise diagonal approximation to data precision matrix
% FORMAT [h] = mda_data_precision (x,mda,load)
%
% 
% h     [P x 1] diagonal of data precision matrix

P = length(x);
dx = 0.001;
h = zeros(P,1);
for p=1:P,
    h(p) = curve(x,mda,load,p,dx);
end

end    

%---------------------------------------------
function [g] = curve (x,mda,load,p,dx)

[x1,x2] = exp_points(x,p,dx);
[L,loglike1,any_nrm] = mda_loglike (x1,mda,load);
[L,loglike2,any_nrm] = mda_loglike (x2,mda,load);

dLdx = (loglike1(:)-loglike2(:))/(2*dx);
g = sum(dLdx.^2);

end

%---------------------------------------------
function [x1,x2] = exp_points(x,p,dx)

x1 = x;
x1(p) = x(p)+dx;
x2 = x;
x2(p) = x2(p)-dx;

end