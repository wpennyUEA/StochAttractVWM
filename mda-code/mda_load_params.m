function [rec,cue] = mda_load_params (mda,params,i)
% Get parameters at specified load
% FORMAT [rec,cue] = mda_load_params (mda,params,i)
% 
% mda       MDA data structure
% params    parameters
% i         load level
%
% rec       parameters of recall model at load level i
% cue       parameters of cue model at load level i

% Load effects for recall attribute
try rec.sigma_e = params.rec.sigma_e(i); catch rec.sigma_e = params.rec.sigma_e(1); end
try rec.sigma = params.rec.sigma(i); catch rec.sigma = params.rec.sigma(1); end
if isfield(params.rec,'sigma_r'),
    try rec.sigma_r = params.rec.sigma_r(i); catch rec.sigma_r = params.rec.sigma_r(1); end
end
if isfield(params.rec,'beta'),  
    try rec.beta = params.rec.beta(i); catch rec.beta = params.rec.beta(1); end
end
if isfield(params.rec,'alpha'),  
    rec.alpha = params.rec.alpha;
end

cue = [];
if mda.opt.mix
    % Load effects for cue attribute
    try cue.sigma_e = params.cue.sigma_e(i); catch cue.sigma_e = params.cue.sigma_e(1); end 
    try cue.sigma = params.cue.sigma(i); catch cue.sigma = params.cue.sigma(1); end
    if isfield(params.cue,'beta'),
        try cue.beta = params.cue.beta(i); catch cue.beta = params.cue.beta(1); end
    end
end
