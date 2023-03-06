function [mda,load] = mda_fit (mda,Y)
% Fit Mixture of Dynamical Attractors Model 
% FORMAT [mda,load] = mda_fit (mda,Y)
%
% mda                   MDA data structure with fields
%   .opt                model options
%
% Y(i)                   data structure at load level i
%
%   .N                   number of items to remember at load level i
%   .T                   [T x 1] vector of delay lengths
%   .S                   [T x N] reported attribute (first column being target)
%   .C                   [T x N] cued attribute (first column being target)
%   .R                   [T x 1] reports (i.e. participant estimates)
%   
% mda                   .params     model parameters
%                       .theta      parameters in latent space
%                       .L          log-likelihood 
%                       .J          log-joint (returned if opt.BayesEst=1)
%                       .Eruns      E of multiple fits 
%                                   (equal to -J for BayesEst, -L otherwise)
%
% load                  data and trial info organised by load
%                       (need this to compute loglike)
%
% Attributes (e.g. location, colour, orientation) must be between 0 and 2*pi

[theta,mda,load] = mda_initialise (mda,Y);
opt = mda.opt;

if opt.BayesEst
    E = mda_neglogjoint (theta,mda,load);
else
    E = mda_negloglike (theta,mda,load);
end
x = theta;

Np = length(theta);

% Initialisations
for i=1:opt.num_restarts,
    theta_init(:,i) = mda_theta_init (opt);
end

if mda.opt.parallel,
    parfor i=1:opt.num_restarts,
        disp(sprintf('Optimisation %d out of %d',i,opt.num_restarts));
        [x_new(:,i),E_new(i),output(i)] = mda_fit_single (theta_init(:,i),mda,load);
    end
else
    for i=1:opt.num_restarts,
        disp(sprintf('Optimisation %d out of %d',i,opt.num_restarts));
        [x_new(:,i),E_new(i),output(i)] = mda_fit_single (theta_init(:,i),mda,load);
    end
end

% Find best fitting model
[E_min,ind] = min(E_new);
x = x_new(:,ind);

R.x = x;
R.E = E_min; 
R.output = output(ind);

mda = mda_fit_report (mda,load,R);
mda.Eruns = E_new;
mda.x_runs = x_new;


