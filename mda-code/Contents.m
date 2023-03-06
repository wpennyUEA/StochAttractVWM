% Matlab Code for [1]
%
% [1] W. Penny (2003) Stochastic Attractor Models of Visual Working Memory
%
% See ../README.m and /demos directory for how to use with your data
%
% Will Penny, School of Psychology, University of East Anglia, UK
% March 2023
%
% All enquiries to: w.penny@uea.ac.uk
%
% -------------------------------------------------------------------------------------
%
%          flow_multi_item - Create flow field for multi-item attractor
%       mda_data_precision - Compute pointwise diagonal approximation to data precision matrix
% mda_data_precision_outer - Compute pointwise outer product approximation to data precision matrix
%            mda_delta_pdf - Delta function PDF at mu
%            mda_diag_hess - Numerical approx to diagonal Hessian using central differences
%                  mda_fit - Fit Mixture of Dynamical Attractors Model
%           mda_fit_report - Complete report of model fitting
%           mda_fit_single - Single fit of MDA model from specified initial condition
%          mda_fit_wrapper - A wrapper for MDA fitting routines that can be called by VB-MFX
%           mda_fp_forward - Forward Solution of Fokker-Planck equation
%              mda_fp_init - Set up Fokker-Planck flow and diffusion matrices
%          mda_fp_matrices - Compute matrices for discretised FP equation
%            mda_group_fit - Fit MDA models to data from a group of subjects using a poor man's MFX
%            mda_group_mfx - Fit MDA models to data from a group of subjects using VB-MFX
%           mda_initialise - Initialise Mixture of Dynamical Attractors Model
%              mda_lat2par - Convert "latent parameters" to parameters
%          mda_load_params - Get parameters at specified load
%              mda_loglike - Compute MDA log-likelihood
%      mda_mixture_density - Compute recall density assuming full MDA model
%          mda_neglogjoint - Compute negative log-joint of MDA model
%           mda_negloglike - Compute negative log-likelihood of MDA model
%              mda_par2lat - Convert parameters to "latent parameters"
%        mda_pop_estimates - Plot distributions of population mean using bootstrapping
%      mda_precompute_flow - Precompute flow function for response-bias model
%         mda_single_trial - Compute recall density for single trial
%              mda_softmax - Softmax function
%           mda_swap_curve - Compute Swap Error Probability for multi-item model
%           mda_theta_init - Randomly initialise model parameters
%             opt_defaults - Set up default parameters for MDA models