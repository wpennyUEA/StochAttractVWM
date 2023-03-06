function [mda] = mda_fit_report (mda,load,R)
% Complete report of model fitting 
% FORMAT [mda] = mda_fit_report (mda,load,R)
%
% INPUTS:
% mda   MDA data structure (from e.g. mda_fit_single)
% load  WM data organised by load
% R     .x          latent params
%       .E          Error 
%       .output     from optimiser
%
% OUTPUTS:
% mda   .aic Akaike Information Criterion
%       .bic Bayesian Information Criterion
%       .F      Laplace approx to model evidence (assuming diag posterior)
%       .Fcorr  Laplace approx to model evidence (assuming full posterior)
%       .Cp  Posterior Covariance Matrix
%       .Rp  Posterior Correlation Matrix
%

x = R.x;
E = R.E;
output = R.output;

Np = length(x);
total_trials = 0;
for i=1:length(load),
    total_trials = total_trials + load(i).Ntrials;
end

% Fitting info (e.g. number of iterations, termination message e.g. solution found)
mda.opt.output = output;
    
% Return estimated parameters
mda.theta = x;
[mda.params,w_lambda,w_alpha,w_w] = mda_lat2par(x,mda);
mda.w_lambda = w_lambda;
mda.w_alpha = w_alpha;
mda.w_w = w_w;

if mda.opt.BayesEst
    mda.L = mda_loglike (mda.theta,mda,load); 
    mda.J = -E;
else
    mda.L = -E;
end
mda.aic = mda.L-Np;
mda.bic = mda.L-0.5*Np*log(total_trials);

if mda.opt.BayesEst
    % Compute posterior covariance and Laplace approx to model evidence
    % See page 217 in Bishop 2006
    disp('Computing Hessian ...');
    
    switch mda.opt.post_precision,
        case 'full',
            % d2L/dtheta_{ij}
            % Computationally expensive and not reliable if local min not found
            A = spm_diff('mda_neglogjoint',x,mda,load,[1 1]);
            mda.Cp = spm_inv(full(A));
            
        case 'diag',
            % Diagonal approximation: ith element (dL/dtheta_i)^2
            h = mda_diag_hess (x,mda,load);
            if any(h<0)
                disp('Error in mda_fit_report:');
                disp('Model has negative covariance');
                keyboard
            end
            mda.Cp = diag(1./h);
            
        case 'pointwise-diag',
            % ith element sum_n (dL_n/dtheta_i)^2
            h = mda_data_precision (x,mda,load);
            % post precision = prior precision + data precision
            h = diag(mda.ipC) + h;
            mda.Cp = diag(1./h);
            
        case 'pointwise-outer',
            % sum_n (dL_n/dtheta)(dL_n/dtheta)^T
            H = mda_data_precision_outer (x,mda,load);
            Lambda = mda.ipC + H;
            mda.Cp = inv(Lambda);
    end
    
    % Account for posterior correlations when approximating model evidence
    mda.Fcorr = mda.J + 0.5*Np*log(2*pi) + 0.5 *spm_logdet(mda.Cp);
    
    % Ignore posterior correlations when approximating model evidence
    mda.F = mda.J + 0.5*Np*log(2*pi) + 0.5 *spm_logdet(diag(diag(mda.Cp)));
    
    % Posterior correlations among parameters in latent space
    sigmap=sqrt(diag(mda.Cp));
    mda.Rp = mda.Cp./(sigmap*sigmap');
    
end
