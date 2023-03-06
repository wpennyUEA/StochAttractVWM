function [opt] = opt_defaults (m,Nloads)
% Set up default parameters for MDA models
% FORMAT [opt] = opt_defaults (m,Nloads)
%
% m         1,2,3 for multi-item, reponse-bias, pure-diffusion resp.
% Nloads    number of different loads (set-sizes)

%------------------------------
% Which type of model?
switch m
    case 1,
        % Multi-Item
        opt.diff_only=0;
        opt.swap=0;
        opt.rb=0;
    case 2,
        % Response-Bias
        opt.diff_only=0;
        opt.swap=1;
        opt.rb=1;
    case 3,
        % Pure-Diffusion
        opt.diff_only=1;
        opt.swap=1;
        opt.rb=0;
end
   
%----------------------------------------------------
% Load effects (2 for an effect of load, 1 otherwise)

opt.Psigma_e = 1;  % encoding
opt.Psigma = Nloads; % delay (diffusion)
opt.Pbeta = 1;      % delay (flow)
opt.Psigma_r = 1;  % recall
opt.load_by_delay=1; % Allow swap and guess rate to vary with load*delay

%--------------------
% Optimisation
opt.num_restarts = 1; % to avoid local maxima
opt.Nbins = 100;
opt.MaxIter = 256;
opt.display = 1;
opt.parallel = 0;
opt.BayesEst = 1;
% opt.post_precision = 'pointwise-diag'; % ith element sum_n (dL_n/dtheta_i)^2
opt.post_precision = 'pointwise-outer'; % sum_n (dL_n/dtheta)(dL_n/dtheta)^T
% opt.post_precision = 'diag'; % ith element (dL/dtheta_i)^2
% opt.post_precision = 'full'; % d2L/dtheta_{ij}


%-----------------------------------------
% Constrain parameter range
opt.lambda_max = 0.2; % max guess rate. Supp Fig 10. in [1] shows max<0.2
opt.alpha_max = 0.2; % max swap rate. Supp Fig 10. in [1] shows max<0.15
opt.beta_max = 15; % flow strength - see Fig 6 in Panichello et al.
opt.sigma_max = 1.5; % diffusion strength
opt.sigma_min = 0.01; % "
opt.min_sd = 0.01; % constraints on encoding and decision noise
opt.max_sd = 0.1; % "

%-----------------
% Various
opt.mix = 0; % multiple attribute model (1) or single attribute (0)
opt.K = 12; % Number of basis functions for RB model
opt.scale_beta = 0; % divide beta by number of items? (1 for yes)