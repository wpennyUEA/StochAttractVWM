function [theta,mda,load] = mda_initialise (mda,Y)
% Initialise Mixture of Dynamical Attractors Model 
% FORMAT [theta,mda,load] = mda_initialise (mda,Y)
%
% INPUTS:
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
%   .R                   [T x 1] location responses 
%
% OUTPUTS:
%
% theta                 (latent) parameters
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

I = length(Y);
mda.I = I;

% Check attributes are between 0 and 2*pi
for i=1:I,
    if sum(any(Y(i).S<0 | Y(i).S>2*pi))
        disp('Error in mda_fit: Y(i).S must be between 0 and 2 pi');
        return
    end
    if isfield(Y(i),'C')
        if sum(any(Y(i).C<0 | Y(i).C>2*pi))
            disp('Error in mda_fit: Y(i).C must be between 0 and 2 pi');
            return
        end
    end
    if any(Y(i).R<0 | Y(i).R>2*pi)
        disp('Error in mda_fit: Y(i).R must be between 0 and 2 pi');
        return
    end
end

opt = mda.opt;
opt.sigma0 = 1.68; % Prior SD
mda.opt = opt;

if opt.load_by_delay
    % Allow swap and guess rate to vary with load*delay
    x=opt.load; y=opt.delay; x0=x-mean(x);y0=y-mean(y);
    I=length(x);T=length(y);
    
    % Create design matrix
    X(:,1) = ones(I*T,1); % average
    X(:,2) = kron(ones(T,1),x0); % effect of load
    X(:,3) = kron(y0,ones(I,1)); % effect of delay
    X(:,4) = X(:,2).*X(:,3); % interaction
    opt.X_load_by_delay = X;
    opt.N_load_by_delay = size(X,2);
    
    % Create indices into output
    r = 1;
    for j=1:T,
        for i=1:I,
            opt.ind_load_by_delay(i,j)=r;
            r = r+1;
        end
    end
    mda.opt = opt;
end

% Initialise parameters
[theta,opt] = mda_theta_init (opt);
mda.opt = opt; % Potentially update prior mean of flow function params for RB models
params = mda_lat2par(theta,mda);

% Record prior mean and Cov
mda.pE = opt.mu0;
mda.pC = opt.sigma0^2*eye(length(mda.pE));
mda.ipC = inv(mda.pC);

% Set up discretisation for Fokker Planck
mda = mda_fp_init (mda);
xc = mda.fp.xc;

% Set up flow fields and initial conditions for each trial
for i=1:I,
    load(i).N = Y(i).N;  % load
    ntrials = length(Y(i).R);
    load(i).Ntrials = ntrials;
    load(i).t = Y(i).T;
    
    for k=1:ntrials,       
        
        load(i).load_index(k)=i; % used for indexing M{i} and A{i,j}
        
        J = length(mda.opt.delay);
        for j=1:J,
            if load(i).t(k) == mda.opt.delay(j)
                load(i).delay_index(k) = j; % used for indexing A{i,j}
            end
        end
        
        c = Y(i).S(k,:); % recall items to remember
        load(i).rec.x(k,:) = c;
        if ~mda.opt.rb,
            % If not response-bias flow, then within-trial flow
            [f,xc] = flow_multi_item (xc,c);
            load(i).rec.f(:,k) = f;
        end
        if opt.mix
            s = Y(i).C(k,:); % cue items to remember
            f = flow_multi_item (xc, s);
            load(i).cue.f(:,k) = f;
            load(i).cue.x(k,:) = s;
            
            [tmp,ind]=min(abs(xc-s(1))); % find bin closest to cue
            load(i).cb(k)=ind;
        end

        [tmp,ind]=min(abs(xc-c(1))); % find bin closest to target
        load(i).target_bin(k)=ind;
        
        r = Y(i).R(k); % participants response
        [tmp,ind]=min(abs(xc-r)); % find bin closest to response
        load(i).rb(k)=ind;
    end
end
