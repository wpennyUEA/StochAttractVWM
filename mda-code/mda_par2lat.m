function [theta] = mda_par2lat (params,mda)
% Convert parameters to "latent parameters" 
% FORMAT [theta] = mda_par2lat (params,mda)
%
% params.rec    recalled attribute:
%               .sigma_e
%               .beta
%               .sigma
%               .sigma_r
%
% params.lambda guess rate (aka non-encoding prob)

opt = mda.opt;

min_sd = opt.min_sd;
max_sd = opt.max_sd;
range_sd = max_sd - min_sd;

if opt.diff_only==0,
    y = params.rec.beta/opt.beta_max;
    theta = logit (y+eps);
end

y = (params.rec.sigma - opt.sigma_min)/(1 - opt.sigma_min);
theta=[theta; logit(y)];

y = (params.rec.sigma_e - min_sd)/range_sd;
theta=[theta; logit(y)];

y = (params.rec.sigma_r - min_sd)/range_sd;
theta=[theta; logit(y)];

if opt.rb | opt.mix
    disp('Error in mda_par2lat.m');
    disp('This function does not work with Response-Bias or Multiple-Attribute models');
end

y = params.lambda/opt.lambda_max;
theta=[theta; logit(y+eps)];

end

%-------------------------------

function [x] = logit(y)

x = log(y./(1-y));

end

