function [mda] = mda_fp_init (mda)
% Set up Fokker-Planck flow and diffusion matrices
% FORMAT [mda] = mda_fp_init (mda)

opt = mda.opt;

% Set up discretisation for Fokker Planck
n = opt.Nbins; % Number of bins (e.g. 100)
[xc,Ddff,Dder] = mda_fp_matrices (n);
mda.fp.xc = xc;
mda.fp.Ddff = Ddff;
mda.fp.Dder = Dder;

if mda.opt.rb
    mda.fp.xb = flow_panichello (mda.fp.xc,opt.K);
end