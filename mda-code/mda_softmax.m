function [y] = mda_softmax(x)
% Softmax function
% FORMAT [y] = mda_softmax(x)

% Keep numbers within machine range
x = x - max(x);

y = exp(x);
y = y/sum(y);