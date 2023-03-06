
clear all
close all

% ----------------------------------------------------------------
% Load your data here
load vwm_data

load = [1 2 4]'; % # different set sizes
delay = [0.5 1 2 4]';  % # different delays in seconds
Nloads = length(load);

% ----------------------------------------------------------------
% Fit Model to Data
%
% See mda_fit.m for required data format 
% and mda_fit.m and mda_fit_report.m for estimated parameters
% and model evidence

model_index = 1; % Multi-Item variant
opt = opt_defaults(model_index,Nloads); 
opt.load = load;
opt.delay = delay;

opt.MaxIter = 16; % change optimisation defaults here e.g. #iterations
    
mda.opt = opt;
tic; mda = mda_fit (mda,Y); toc

% ----------------------------------------------------------------
% Parameter Estimates

disp(' ');
disp('Loads:');
disp(load);

disp('Delays (secs):');
disp(delay);

disp('Attractor Strength:');
disp(mda.params.rec.beta);

disp('Diffusion Strength at different Loads:');
disp(mda.params.rec.sigma);

disp('Guess Rates at Load x Delay levels:');
disp(mda.params.lambda);