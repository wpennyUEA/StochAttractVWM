function [group,hier,D,F] = mda_group_mfx (mda_org,sub,ind)
% Fit MDA models to data from a group of subjects using VB-MFX
% FORMAT [group,hier,D,F] = mda_group_mfx (mda_org,sub,ind)
%
% INPUTS:
% mda_org       MDA data structure
% sub           sub(i).Y is data for subject i
% ind           indices of parameters to "average" e.g. [6,7,8,9]
%               IGNORED FOR NOW
%
% OUTPUTS:
% group         group.mda{i} is fitted MDA model for subject i
% hier          2nd level model
% D             first level estimates with and without MFX
% F             overall log model evidence

%----------------------------------------------------------------
% Set up group-level prior precisions

[theta,mda] = mda_initialise (mda_org,sub(1).Y);
Np = length(theta);
all = [1:Np];

precision='vague';
switch precision
    case 'within'
        % same precision as within-subject
        mean_prec = 1./diag(mda.pC);
        var_prec = 0.01*ones(Np,1);
        hier.b0 = mean_prec./var_prec;
        hier.a0 = hier.b0.*mean_prec;
        
    case 'low-and-high',
        % low-precision settings
        mean_prec=1;
        var_prec=10;
        hier.b0(all) = mean_prec/var_prec;
        hier.a0(all) = hier.b0(all)*mean_prec;
        
        % high-precision settings e.g. for RB flow-function params
        % i.e. not much variability about group mean
        mean_prec=10;
        var_prec=100;
        hier.b0(ind) = mean_prec/var_prec;
        hier.a0(ind) = hier.b0(ind)*mean_prec;
        
    case 'vague',
        mean_prec = 10*ones(Np,1);
        var_prec = 1000*ones(Np,1);
        hier.b0 = mean_prec./var_prec;
        hier.a0 = hier.b0.*mean_prec;
end


%------------------------------------------------------------------
% Initialise models
Nsub = length(sub);
for i=1:Nsub,
    [theta(:,i),mda,load] = mda_initialise (mda_org,sub(i).Y);
    
    % Create info for mda_fit_wrapper
    model{i}.mda = mda;
    model{i}.load = load;
    model{i}.Y = sub(i).Y;
    
    model{i}.mu = mda.pE;
    model{i}.S = inv(mda.pC);
end

opt =[]; % accept VB-MFX defaults
opt.max_subj_its = 240;
%opt.max_subj_its = 120;
%opt.max_subj_its = 16;
[model,hier,D,F] = vbmfx (model,'mda_fit_wrapper',hier,opt);

for i=1:Nsub,
    group.mda{i} = model{i}.mda;
end
