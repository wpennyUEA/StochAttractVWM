function [L,loglike,any_nrm] = mda_loglike (theta,mda,load)
% Compute MDA log-likelihood
% FORMAT [L,loglike,any_nrm] = mda_loglike (theta,mda,load)
%
% theta     latent parameters
% mda       MDA data structure
% load      
%           load(i).t(k) delay length
%           load(i).rb(k) participant's response bin
%           load(i).cb(k) participant's cue bin
%           load(i).load_index(k) (=i) used to index M{i} and A{i,j}
%           load(i).delay_index(k) (=j) used to index A{i,j}
%
%           For recalled attribute:
%               load(i).rec.f(:,k) flow field for kth trial at load level i
%               load(i).rec.x(k,:) vector of recall-attribute items
%
%           For cued attribute:
%               load(i).cue.f(:,k) flow field for kth trial at load level i
%               load(i).cue.x(k,:) vector of cue-attribute items
%
% L         log-likelihood
% loglike   loglike(i,k) is log likelihood for trial k at load level i
% any_nrm   1 if any negative entries removed in P (see mda_single_trial), 0 otherwise
%
% For response-bias or pure-diffusion (f=0) models 
%   M{i}=sigma(i)*D_n -beta(i)*Da*diag(f) and
%   A{i,j} = exp(M{i} t(j)) 
% are condition-specific (rather than trial-specific) so can be pre-computed 
% to save computation

tau=0; % delay length
params = mda_lat2par (theta,mda);

mda.fp.Aenc = expm(mda.fp.Ddff * params.rec.sigma_e);
mda.fp.Adec = expm(mda.fp.Ddff * params.rec.sigma_r);
            
if mda.opt.rb | mda.opt.diff_only
    if mda.opt.rb
        % Flow function is constant across trials for response-bias model
        df = mda.fp.xb*params.rec.w;
        df = df / max(abs(df));
    end
    % Pre-compute transition matrix for each combination of load and delay
    I = length(mda.opt.load);
    J = length(mda.opt.delay);
    
    for i=1:I,
        if mda.opt.rb
            if length(params.rec.beta) > 1
                beta = params.rec.beta(i);
            else
                beta = params.rec.beta;
            end
        end
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
            A1 = mda.fp.Aenc;
            A2 = expm(M*mda.opt.delay(j));
            A3 = mda.fp.Adec;
            mda.fp.A{i,j} = A1*A2*A3;
            
            %mda.fp.A{i,j} = expm(M*mda.opt.delay(j));
        end
    end
end

any_nrm=0;

% Uniform density for non-encoding trials
Pnon = ones(mda.opt.Nbins,1)/mda.opt.Nbins;

% Loop over loads
L=0;
for i=1:mda.I,
    ntrials=length(load(i).rb);
    [rec,cue] = mda_load_params (mda,params,i);
    
    for k=1:ntrials,
        trial.num_items = load(i).N;
        trial.t = load(i).t(k);
        trial.i = load(i).load_index(k);
        trial.j = load(i).delay_index(k);
        rec.x = load(i).rec.x(k,:);
        
        if mda.opt.rb
            trial.df = df;
            trial.A = mda.fp.A {trial.i,trial.j};
        else
            trial.df = load(i).rec.f(:,k);
        end
        
        if mda.opt.mix
            trial.df_cue = load(i).cue.f(:,k);
            cue.x = load(i).cue.x(k,:);
            cue.bin = load(i).cb(k);
            [P,nrm] = mda_single_trial (mda,trial,rec,cue,params.lambda);
        else
            [P,nrm] = mda_single_trial (mda,trial,rec);
            
            % Account for possibility item is not encoded
             if mda.opt.load_by_delay
                lambda = params.lambda(trial.i,trial.j);
             else
                lambda = params.lambda;
             end
             P = lambda*Pnon+(1-lambda)*P;
        end
        
        if nrm, any_nrm=1; end
        
        response_bin = load(i).rb(k);
        loglike(i,k) = log(P(response_bin));
        if isnan(loglike(i,k)) | ~isreal(loglike(i,k))
            keyboard
        end
    end
end
L = sum(sum(loglike));