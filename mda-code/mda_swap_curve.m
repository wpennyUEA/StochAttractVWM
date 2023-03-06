function [se] = mda_swap_curve (x,xc,g,P)
% Compute Swap Error Probability for multi-item model
% FORMAT [se] = mda_swap_curve (x,xc,g,P)
%
% x     recall attribute items, first entry is target, others non-targets
% xc    discrete space of possible attribute values
% g     flow function
% P     probability density over xc
%
% se    swap error
%
% Unstable Fixed Points assumed half way between stable ones (as in MI
% model)

N = length(x); % Number of items to maintain

if N == 1,
    se=0;
else
    % Find unstable fixed points (UFPs)
    [sx,ind] = sort(x);
    for i=1:N-1,
        mid(i) = mean([sx(i),sx(i+1)]);   
    end
    delta = 0.5*(sx(1)+2*pi-sx(N));
    mid(N) = mod(sx(N) + delta,2*pi);
    
    % Check identification of UFPs
%     figure;
%     plot(xc,g);
%     hold on
%     for i=1:N,
%         plot(sx(i),0,'ro');
%         plot(mid(i),0,'rx');
%     end
%     grid on
    
    % Compute basin of attraction for ith item
    for i=1:N,
        if i==1
            lower(i) = mid(N);
        else
            lower(i) = mid(i-1);
        end
        upper(i) = mid(i);
    end
    
    % Remove target
    ind=find(sx==x(1));
    lower(ind)=[]; upper(ind)=[];
    
    % Compute total integral between upper and lower bounds
    se = 0;
    for i=1:N-1,
        [tmp,low_ind]=min(abs(xc-lower(i)));
        [tmp,high_ind]=min(abs(xc-upper(i)));
        
        if upper(i) > lower(i)
            se = se + sum(P(low_ind:high_ind));
        else
            se = se + sum(P(low_ind:end));
            se = se + sum(P(1:high_ind));
        end
    end
end
        
        