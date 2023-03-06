function [H] = mda_data_precision_outer (x,mda,load)
% Compute pointwise outer product approximation to data precision matrix
% FORMAT [H] = mda_data_precision_outer (x,mda,load)
%
% 
% h     [P x 1] data precision vector

P = length(x);
dx = 0.001;
for p=1:P,
    dLdx(:,p) = grad(x,mda,load,p,dx);
end

N = size(dLdx,1);
H = zeros(P,P);
for n=1:N,
    H = H + dLdx(n,:)'*dLdx(n,:);
end

end    

%---------------------------------------------
function [dLdx] = grad (x,mda,load,p,dx)

[x1,x2] = exp_points(x,p,dx);
[L,loglike1,any_nrm] = mda_loglike (x1,mda,load);
[L,loglike2,any_nrm] = mda_loglike (x2,mda,load);

dLdx = (loglike1(:)-loglike2(:))/(2*dx);

end

%---------------------------------------------
function [x1,x2] = exp_points(x,p,dx)

x1 = x;
x1(p) = x(p)+dx;
x2 = x;
x2(p) = x2(p)-dx;

end