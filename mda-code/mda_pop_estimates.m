function [] = mda_pop_estimates (x,names)
% Plot distributions of population mean using bootstrapping
% FORMAT [] = mda_pop_estimates (x,names)
%
% x         [D x N] matrix where D is number of variabes, N subjects
% names     {d} name of dth variable, d = 1..D

x = x';
[N,D] = size(x);
%if D>N, x=x'; [N,D] = size(x); end

B = 1000;
% bootstrap resampling indices
ind=ceil(rand(N,B)*N);
for b=1:B,
    
end
for d=1:D,
    a = x(:,d);
    for b=1:B,
        S(b,d) = mean(a(ind(:,b)));
    end
end

figure
%violinplot(S,names);
violinplot(S,names,'ShowData',logical(0));
%boxplot(S);
grid on

if D==2
    diff=S(:,2)-S(:,1);
    pr=length(find(diff>0))/B;
    pr = max([pr,1-pr]);
    disp(sprintf('Proportion of samples passing sign test = %1.3f',pr));
end