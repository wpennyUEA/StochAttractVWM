function [mda,df] = mda_precompute_flow (mda,params)
% Precompute flow function for response-bias model
% FORMAT [mda,df] = mda_precompute_flow (mda,params)

% Flow function is constant across trials for response-bias model
df = mda.fp.xb*params.rec.w;
df = df / max(abs(df));

% Pre-compute transition matrix for each combination of load and delay
I = length(mda.opt.load);
J = length(mda.opt.delay);

for i=1:I,
    beta = params.rec.beta(i);
    sigma = params.rec.sigma(i);
    num_items = mda.opt.load(i);
    if mda.opt.diff_only == 1,
        % Diffusion term only
        M = mda.fp.Ddff * sigma;
    else
        if mda.opt.scale_beta, beta_s = beta/num_items; else beta_s=beta; end
        M = mda.fp.Ddff * sigma  - mda.fp.Dder * diag(df) * beta_s;
    end
    
    for j=1:J,
        %*************************************************
        % Assumes sigma_e and sigma_r don't vary with load
        % ************************************************
        A1 = expm(mda.fp.Ddff * params.rec.sigma_e);
        A2 = expm(M*mda.opt.delay(j));
        A3 = expm(mda.fp.Ddff * params.rec.sigma_r);
        mda.fp.A{i,j} = A1*A2*A3;
        
        %mda.fp.A{i,j} = expm(M*mda.opt.delay(j));
    end
end