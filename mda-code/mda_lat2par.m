function [params,w_lambda,w_alpha,w_w] = mda_lat2par (theta,mda)
% Convert "latent parameters" to parameters
% FORMAT [params,w_lambda,w_alpha,w_w] = mda_lat2par (theta,mda)
%
% params        .lambda     guess rate
%
% params.rec    recalled attribute:
%               .sigma_e
%               .beta
%               .sigma
%               .sigma_r
%               .w
%               .alpha      swap rate
%
% params.cue    cued attribute:
%               .sigma_e
%               .beta
%               .sigma
%
%               Use these to test for main effects and interactions:
% w_lambda      latent params for guess rate effects 
% w_alpha       latent params for swap rate effects
% w_w           latent params for flow function

w_lambda = [];
w_alpha = [];
w_w = [];

opt = mda.opt;

min_sd = opt.min_sd;
max_sd = opt.max_sd;
range_sd = max_sd - min_sd;

if opt.diff_only==0,
    params.rec.beta = sigmoid(theta(1:opt.Pbeta))*opt.beta_max;
    theta(1:opt.Pbeta)=[];
end

params.rec.sigma = opt.sigma_min + sigmoid(theta(1:opt.Psigma))*(opt.sigma_max - opt.sigma_min);
theta(1:opt.Psigma)=[];

params.rec.sigma_e = min_sd+sigmoid(theta(1:opt.Psigma_e))*range_sd;
theta(1:opt.Psigma_e)=[];

params.rec.sigma_r = min_sd+sigmoid(theta(1:opt.Psigma_r))*range_sd;
theta(1:opt.Psigma_r)=[];

if opt.rb
    params.rec.w = mda_softmax(theta(1:opt.K));
    w_w = theta(1:opt.K);
    %params.rec.w = theta(1:opt.K);
    theta(1:opt.K)=[];
end

if opt.mix==0, 
    if opt.load_by_delay
        % Guess rate
        w_lambda = theta(1:opt.N_load_by_delay);
        a_lambda = opt.X_load_by_delay*w_lambda;
        theta(1:opt.N_load_by_delay)=[];
        
        if opt.swap
            % Swap rate
            w_alpha = theta(1:opt.N_load_by_delay);
            a_alpha = opt.X_load_by_delay*w_alpha;
            theta(1:opt.N_load_by_delay)=[];
        end
        
        for i=1:length(opt.load),
            for j=1:length(opt.delay),
                r = opt.ind_load_by_delay(i,j);
                params.lambda(i,j) = sigmoid (a_lambda(r))*opt.lambda_max;
                if opt.swap
                    params.rec.alpha(i,j) = sigmoid(a_alpha(r))*opt.alpha_max;
                end
            end
        end
    else
        params.lambda = sigmoid(theta(1))*opt.lambda_max;
        theta(1) = [];
        if opt.swap
            params.rec.alpha = sigmoid(theta)*opt.alpha_max;
        end
    end
    
    return; 
end

% ---------------------------------------------------
% For MDA models - add parameters for cued attribute

if opt.diff_only==0,
    params.cue.beta = sigmoid(theta(1:opt.Pbeta));
    theta(1:opt.Pbeta)=[];
end
params.cue.sigma = sigmoid(theta(1:opt.Psigma));
theta(1:opt.Psigma)=[];

params.cue.sigma_e = min_sd+sigmoid(theta(1:opt.Psigma_e))*range_sd;
theta(1:opt.Psigma_e)=[];
    
params.lambda = sigmoid(theta)*opt.lambda_max;


end

%-------------------------------

function [y] = sigmoid(x)

y = 1./(1+exp(-x));

end

