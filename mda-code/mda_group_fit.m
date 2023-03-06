function [group] = mda_group_fit (mda_org,sub,ind)
% Fit MDA models to data from a group of subjects using a poor man's MFX
% FORMAT [group] = mda_group_fit (mda_org,sub,ind)
%
% mda_org       MDA data structure
% sub           sub(i).Y is data for subject i
% ind           indices of parameters to average e.g. [6,7,8,9]
%
% group         group.mda{i} is fitted MDA model for subject i

initial_its = 16;
subs_its = 16;
group_its = 4;

mda_org.opt.MaxIter = initial_its;

% Rough individual subject fits
Nsub = length(sub);
for i=1:Nsub,
    [mda{i},load{i}] = mda_fit (mda_org,sub(i).Y);
    theta(:,i) = mda{i}.theta;
end

Np = length(mda{1}.theta);
pE = zeros(Np,1);

for it = 1:group_its,
    
     disp(sprintf('Group Optimisation Iteration %d ...',it));
     
     % Compute mean and SEM of (selected) parameter estimates over group 
     for i=1:Nsub,
         x(:,i) = theta(ind,i);
     end
     m = mean(x')';
     s = std(x')'/Nsub;
     v = s.^2;
     
     for i=1:Nsub,
         
        % Change initialisation
        theta_init(:,i) = theta(:,i);
        theta_init(ind,i) = m;
        
        % Change prior mean
        mda{i}.pE (ind) = m;
        
        % Change prior cov
%         vold = diag(mda{i}.pC);
%         vnew = vold;
%         vnew(ind) = v;
%         mda{i}.pC = diag(vnew);
     end
     
     % Refit models with group priors and new initialisation
     for i=1:Nsub,
         mda{i}.opt.MaxIter = subs_its;
         [theta(:,i),E_new(i),output(i)] = mda_fit_single (theta_init(:,i),mda{i},load{i});
        %mda{i} = mda_fit (mda{i},sub(i).Y);
     end
     
end

% Complete model fit reports e.g. compute evidence, format params etc
for i=1:Nsub,
    R.x = theta(:,i);
    R.E = E_new(i);
    R.output = output(i);
    mda{i} = mda_fit_report (mda{i},load{i},R);
end

group.mda = mda;