function [P] = mda_delta_pdf (mda,mu)
% Delta function PDF at mu
% FORMAT [P] = mda_delta_pdf (mda,mu)

x = mda.fp.xc;
[tmp,ind]=min(abs(x-mu)); % find bin closest to x0
n = length(x);
P = zeros(n,1);
P(ind) = 1;