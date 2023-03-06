function [E] = mda_neglogjoint (theta,mda,load)
% Compute negative log-joint of MDA model
% FORMAT [E] = mda_neglogjoint (theta,mda,load)

pE = mda.pE;
pC = mda.pC;
ipC = diag(1./diag(pC));

% Log Likelihood
L1 = mda_loglike (theta,mda,load); 

% Log Prior
Np = length(theta);
e = theta-pE;
L0 = -0.5*e'*ipC*e;
L0 = L0 - 0.5*sum(log(diag(pC))) -0.5*Np*log(2*pi);

% Old code:
% mu = mda.opt.mu0;
% s0 = mda.opt.sigma0;
% L0 = -sum((theta-mu).^2)/(2*s0^2);
% L0 = L0 -Np*log(s0)-0.5*Np*log(2*pi);

E = -L1-L0;

