function [f,xc] = flow_multi_item (xc,c)
% Create flow field for multi-item attractor
% FORMAT [f,xc] = flow_multi_item (xc,c)
%
% xc    [Nbins x 1] bin locations
% c     [N x 1] items to remember (locations of stable fixed points)
%
% f     [Nbins x 1] flow vector
%
% Adapt phase velocity to ensure single cycle between N fixed points

N = length(c);

[cs,ind] = sort(c);
for j=1:N,
    % Find bin nearest cs(j)
    [tmp,ind]=min(abs(xc-cs(j)));
    bin(j)=ind;
end

n = length(xc);
dx = xc(2)-xc(1);

f = zeros(n,1);
for j=1:N-1,
    start = bin(j);
    stop = bin(j+1);
    i = [start:stop];
    xstart = xc(start);
    xstop = xc(stop);
    v = (2*pi)/(xstop-xstart); % speed
    x = [xstart:dx:xstop]; % domain
    phi = v*(x-xstart);
    f(i) = -sin(phi);
end

% Work out velocity for final leg
xstart = xc(bin(N));
xstop = xc(bin(1))+2*pi;
v = (2*pi)/(xstop-xstart);
    
% From bin(N) to end
i = [bin(N):n];
x = xc(i);
phi = v*(x-x(1));
f(i) = -sin(phi);
x0 = x(1);

% From 1 to bin(1)
phi_old=phi(end);
i = [1:bin(1)];
x = xc(i)+2*pi;
phi = v*(x-x0);
f(i) = -sin(phi);


