function [E] = mda_negloglike (theta,mda,load)
% Compute negative log-likelihood of MDA model
% FORMAT [E] = mda_negloglike (theta,mda,load)

L = mda_loglike (theta,mda,load); 
E = -L;

