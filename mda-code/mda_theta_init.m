function [theta,opt] = mda_theta_init (opt)
% Randomly initialise model parameters
% FORMAT [theta,opt] = mda_theta_init (opt)

if opt.diff_only==0,
    theta = opt.sigma0*randn(opt.Pbeta,1);
    mu = zeros(opt.Pbeta,1);
else
    theta=[];
    mu = [];
end

theta = [theta; opt.sigma0*randn(opt.Psigma,1)]; 
mu = [mu; zeros(opt.Psigma,1)];
theta = [theta; opt.sigma0*randn(opt.Psigma_e,1)];
mu = [mu; zeros(opt.Psigma_e,1)];
theta = [theta; opt.sigma0*randn(opt.Psigma_r,1)];
mu = [mu; zeros(opt.Psigma_r,1)];

if opt.rb
    if isfield(opt,'m_w')
        theta = [theta; opt.m_w];
        mu = [mu; opt.m_w];
    else
        theta = [theta; opt.sigma0*randn(opt.K,1)];
        mu = [mu; zeros(opt.K,1)];
    end
end

if opt.mix==0, 
    if opt.load_by_delay
        % Allow guess rate (and maybe swap rate) to vary with load*delay
        Ntheta = opt.N_load_by_delay;
        if opt.swap
            Ntheta = Ntheta*2; % Twice as many params if we're modelling swaps
        end
        theta = [theta; opt.sigma0*randn(Ntheta,1)];
        mu = [mu; zeros(Ntheta,1)];
    else
        % for non-encoding prob
        theta = [theta; opt.sigma0*randn(1,1)];
        mu = [mu; 0];
        if opt.swap
            theta = [theta; opt.sigma0*randn(1,1)];
            mu = [mu; 0];
        end
    end
    
    opt.mu0 = mu;
    return; 
end

% ---------------------------------------------------
% For MDA models - add parameters for cued attribute

if opt.diff_only==0,
    theta = [theta; opt.sigma0*randn(opt.Pbeta,1)];
    mu = [mu; zeros(opt.Pbeta,1)];
end
theta = [theta; opt.sigma0*randn(opt.Psigma,1)];
mu = [mu; zeros(opt.Psigma,1)];
theta = [theta; opt.sigma0*randn(opt.Psigma_e,1)];
mu = [mu; zeros(opt.Psigma,1)];

theta = [theta; opt.sigma0*randn(1,1)]; % for prob of encoding failure, pi
mu = [mu; 0];

opt.mu0 = mu;
