function [xc,Ddff,Dder] = mda_fp_matrices (n)
% Compute matrices for discretised FP equation
% FORMAT [xc,Ddff,Dder] = mda_fp_matrices (n)
%
% n         number of bins
%
% xc        [n x 1] bin locations
% Ddff      [n x n] diffusion matrix
% Dder      [n x n] flow matrix

dx = 2*pi/n; %bin width
%xe = linspace(-pi,pi,n+1); %bin edges
xe = linspace(0,2*pi,n+1); %bin edges
xc = (xe(1:n) + dx/2)'; %bin centers

cDiff = 1/(dx^2*2);  % scale factor for diffusion
cDrft = 1/(dx*2);   % scale factor for drift

% diffusion matrix
Ddff = [-2 * cDiff, cDiff, zeros(1,n-3), cDiff];
Ddff = toeplitz([Ddff(1) fliplr(Ddff(2:end))],Ddff);

% drift matrix
Dder = [0, cDrft, zeros(1,n-3), -cDrft];
Dder = toeplitz([Dder(1) fliplr(Dder(2:end))],Dder);


