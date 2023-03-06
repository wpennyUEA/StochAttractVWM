function [P,nrm] = mda_fp_forward (M,t,P0,A)
% Forward Solution of Fokker-Planck equation
% FORMAT [P,nrm] = mda_fp_forward (M,t,P0,A)
%
% M,t,P0    
% A         pre-computed expm(M*t) 
%
% P         expm(M * t) * P0
% nrm       1 indicates negative values removed 
%           (if so, they are replaced by 0, then P renormalised)

if nargin < 4
    P = expm(M * t) * P0;
else
    P = A * P0;
end

if any(P<0)
    ind = find(P<0);
    P(ind) = 0;
    P = P/sum(P);
    nrm = 1;
else
    nrm = 0;
end

