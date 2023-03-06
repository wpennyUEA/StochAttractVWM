function [model] = mda_fit_wrapper (model,opt)
% A wrapper for MDA fitting routines that can be called by VB-MFX
% FORMAT [model] = mda_fit_wrapper (model,opt)
%
% INPUT:
%
% model.mu       prior mean  (v in paper)
%      .S        prior precision  (lambda in paper)
%      .Y        data
%      .mda      MDA data structure
%      .load     data organised by load and other precomputed quantities
%
% opt           optimisation options (see vbmfx.m)
%               .maxits
%
% OUTPUT:
%
% model.w        posterior mean 
% model.R        posterior precision 
% model.F        approximation to model evidence
% model.mda      MDA data structure


mda = model.mda;
load = model.load;

mda.pE = model.mu;
mda.pC = inv(model.S);
mda.ipC = model.S;

%theta_init = mda.pE;
theta_init = spm_normrnd(mda.pE, mda.pC,1);

mda.opt.MaxIter = opt.maxits;
[theta,E_new,output] = mda_fit_single (theta_init,mda,load);

A.x = theta;
A.E = E_new;
A.output = output;
mda = mda_fit_report (mda,load,A);

model.w = mda.theta;
model.R = inv(mda.Cp);
model.F = mda.F;
model.mda = mda;
