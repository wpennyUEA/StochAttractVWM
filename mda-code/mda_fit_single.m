function [x,E,output] = mda_fit_single (theta,mda,load)
% Single fit of MDA model from specified initial condition
% FORMAT [x,E,output] = mda_fit_single (theta,mda,load)
%
% theta         initial condition
% mda           model
% load          data
%
% x             estimated params
% E             value of error function

% The default optimisation 'fminsearch' uses the Nelder-Mead Simplex method
% The 'fminunc' algorithm uses a Trust Region method if you supply the
% gradient but, if not, will use a Quasi-Newton method.

try display = mda.opt.display; catch display=0; end

Np=length(theta);
options.MaxFunEvals=128;
if display
    options = optimset(options,'Display','iter');
else
    options = optimset(options);
end

if mda.opt.BayesEst
    costfun = 'mda_neglogjoint';
else
    costfun = 'mda_negloglike';
end

%opt = 'fminsearch';
algo = 'fminunc';
switch algo
    case 'fminunc',
%         options.TolFun=0.01; % Reduce this to get a better solution
%         options.TolX =1e-4;
        options.TolFun=1e-4; 
        options.TolX =1e-6;
        options.MaxIter = mda.opt.MaxIter;
        options.MaxFunEvals=512*Np;
        
        
        % As used in Panichello ...
%         options.TolX =1e-8;
%         options.MaxFunEvals=10000;
        
        [x,E,eflag,output]=fminunc(costfun,theta,options,mda,load);
        
    case 'fminsearch',
        options.TolFun=0.1;
        options.TolX =1;
        [x,E]=fminsearch(costfun,theta,options,mda,load);
    otherwise
        disp('Error: Unknown optimiser in mda_fit.m');
        return
end