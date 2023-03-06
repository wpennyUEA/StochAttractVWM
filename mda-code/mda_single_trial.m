function [P,nrm] = mda_single_trial (mda,trial,rec,cue,pi0)
% Compute recall density for single trial
% FORMAT [P,nrm] = mda_single_trial (mda,trial,rec,cue,pi0)
% 
% mda       MDA data structure
%               .opt.Nbins  Number of bins for Fokker-Planck discretisation
%               .opt.diff_only 1 for Pure-Diffusion
%               .opt.mix    1 for Multiple Attributes Model
%               .opt.swap   1 for Additive swap error term
%
% trial     trial specific info
%               .df         flow vector
%               .t          delay length(s)
%               .num_items  number of items
%               .i          load index 
%               .j          delay index
%
% rec       parameters of recalled attribute
%               .x          recall attribute vector
%               .beta       flow strength
%               .sigma      diffusion strength
%               .alpha      swap error prob
%               .sigma_e    SD of encoding noise
%               .sigma_r    SD of decision noise
%
% cue       parameters of cued attribute
%               .x      cue attribute vector
%
% pi0       non-encoding probability
%
% P         [Nbins x 1] recall density
% nrm       1 if any negative entries removed, 0 otherwise

Ddff = mda.fp.Ddff;
Dder = mda.fp.Dder;

df = trial.df;
t = trial.t;
num_items = trial.num_items;

if ~mda.opt.rb
    % Markov matrix
    if mda.opt.diff_only==1,
        % Diffusion term only
        Mrec = Ddff * rec.sigma;
    else
        if mda.opt.scale_beta, beta_s = rec.beta/num_items; else beta_s=rec.beta; end
        Mrec = Ddff * rec.sigma  - Dder * diag(df) * beta_s;
    end
    % Note: we optionally divide by num_items to correct for the increase
    % in the gradient of df (and therefore strength of attractor) with number of items
end

if ~mda.opt.mix
    % Single Dynamic Attractor (SDA)
    
    mu = rec.x(1);
    P0 = mda_delta_pdf (mda,mu); 
    m = length(rec.x(2:end));
    
    if m > 1 & mda.opt.swap
        % Additional swap error term
        if mda.opt.load_by_delay
            alpha = rec.alpha(trial.i,trial.j);
        else
            alpha = rec.alpha;
        end
        P0 = (1-alpha)*P0;
        for i=1:m,
            mu = rec.x(i+1);
            P0 = P0 + (alpha/m)*mda_delta_pdf (mda,mu);
        end
    end
    
    % Density after delay t
    if mda.opt.rb | mda.opt.diff_only
        % Use pre-computed transition matrix
        A = mda.fp.A{trial.i,trial.j};
        [P,nrm] = mda_fp_forward ([],[],P0,A);
    else
        if isfield(mda.fp,'Aenc')
            % Use Precomputed Flow Fields for Encoding and Decision
            
            % Add encoding noise
            Aenc = mda.fp.Aenc;
            [P0,nrm] = mda_fp_forward ([],[],P0,Aenc);
            
            % Maintenance
            [P,nrm] = mda_fp_forward (Mrec,t,P0);
            
            % Add decision noise
            Adec = mda.fp.Adec;
            [P,nrm] = mda_fp_forward ([],[],P,Adec);
        else
            % Add encoding noise
            R = Ddff * rec.sigma_e;
            [P0,nrm] = mda_fp_forward (R,1,P0);
            
            % Maintenance
            [P,nrm] = mda_fp_forward (Mrec,t,P0);
            
            % Add decision noise
            R = Ddff * rec.sigma_r;
            [P,nrm] = mda_fp_forward (R,1,P);
        end
        
        
    end
    
else
    % Mixture of Dynamical Attractors (MDA)
    
    if mda.opt.diff_only==1,
        % Diffusion term only
        Mcue = Ddff * cue.sigma;
    else
        if mda.opt.scale_beta, beta_s = cue.beta/num_items; else beta_s=cue.beta; end
        Mcue = Ddff * cue.sigma  - Dder * diag(trial.df_cue) * beta_s;
        keyboard
    end
    
    % Uniform densities for non-encoding slot
    Nbins = length(mda.fp.xc);
    Prec(:,1) = ones(Nbins,1)/Nbins;
    cue_like(1) = 1/Nbins;
    
    % Priors
    pr(1) = pi0;
    pr(2:num_items+1) = (1-pr(1))*ones(num_items,1)/num_items;
    
    % get item-specific recall densities
    for j=1:num_items,
        m = rec.x(j);
        P0 = mda_delta_pdf (mda,m); 
        
        % Add encoding noise
        R = Ddff * rec.sigma_e;
        [P0,nrm] = mda_fp_forward (R,1,P0);
    
        % Density after delay t
        %Prec(:,j+1) = expm(Mrec * t) * P0;
        [Prec(:,j+1),nrm] = mda_fp_forward (Mrec,t,P0);
    end
    
    % get item-specific cue densities
    for j=1:num_items,
        
        m = cue.x(j);
        P0 = mda_delta_pdf (mda,m); 
        
        % Add encoding noise
        R = Ddff * cue.sigma_e;
        [P0,nrm] = mda_fp_forward (R,1,P0);
        
        % Density after delay t
        % Pcue = expm(Mcue * t) * P0;
        [Pcue,nrm] = mda_fp_forward (Mcue,t,P0);
    
        if any(Pcue<0)
            % If sigma is small in relation to (2*pi)/Nbins we may get
            % non-zero entries in P. If so, remove them.
            ind = find(Pcue<0);
            Pcue(ind) = 0;
            Pcue = Pcue/sum(Pcue);
        end
    
        cue_like(j+1) = Pcue (cue.bin);
    end
    
    % get posterior density over slots given cue
    slot_post = cue_like.*pr;
    slot_post = slot_post/sum(slot_post);
    
    % get recall density given cue
    P = Prec*slot_post(:);
    
    % Add decision noise
    R = Ddff * rec.sigma_r;
    [P,nrm] = mda_fp_forward (R,1,P);
    

end


