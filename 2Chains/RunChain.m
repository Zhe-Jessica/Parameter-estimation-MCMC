close all; clear all; clc;
%% Define constants
n = [314.755, 726.0529, 307.3519];    % beta distribution 
x0 = [90000, 0.5, 0.5];                                                                           % initial values of parameters 

%% MCMC preparation
nIter = 5000;                                                      % iterations of the Metropolis-Hastings algorithm
% PI, GOR, WC
record = zeros(nIter,3);
reject = zeros(nIter,3);
record(1,:) = x0';
[P_current,y_current] = posterior(x0(1), x0(2), x0(3));
Results(:,:,1) = y_current;
RejectRate = 0;

%% MCMC
% Core of the algorithm
for k = 2:nIter
    k
    current = record(k - 1,:);
    % transfer to (0,1)
    ValueToAlpha = (current(1)-10^4)/(10^5-10^4);
    transfer = betarnd(ValueToAlpha*n(1),(1-ValueToAlpha)*n(1));
    proposed_PI = (10^5-10^4)*transfer+10^4;

    proposed_GOR = betarnd(current(2)*n(2),(1-current(2))*n(2));
    proposed_WC = betarnd(current(3)*n(2),(1-current(3))*n(2));

    % Metropolis-Hustings evaluation
    [P_proposed,y_proposed] = posterior(proposed_PI, proposed_GOR, proposed_WC);
    A = P_proposed -  P_current + log(betapdf(ValueToAlpha,n(1)*transfer,n(1)*(1-transfer))) +log(betapdf(current(2),n(2)*proposed_GOR,n(2)*(1-proposed_GOR))) - log(betapdf(proposed_GOR,n(2)*current(2),n(2)*(1-current(2))))+log(betapdf(current(3),n(3)*proposed_WC,n(3)*(1-proposed_WC))) - log(betapdf(proposed_WC,n(3)*current(3),n(3)*(1-current(3))));

    if(log(rand) < A)                                                  %accept
        if transfer == 0 | transfer == 1 ;
            proposed_PI = current(1);
        end
        if proposed_GOR == 0 | proposed_GOR == 1 ;
            proposed_GOR = current(2);
        end
        if proposed_WC == 0|proposed_WC == 1;
            proposed_WC = current(3);
        end
         record(k,:) = [proposed_PI; proposed_GOR; proposed_WC];
         P_current = P_proposed;
         y_current = y_proposed;
    else                                                                        %reject
         record(k,:) = current';
         reject(k,:) = [proposed_PI; proposed_GOR; proposed_WC];
         RejectRate = RejectRate +1;
    end
    Results(:,:,k) = y_current;                     % save ye

end

save("Chain1.mat")

figure
tiledlayout(3,2)
nexttile
histogram(record(:,1))
title("PI0 = "+x0(1)+"")
nexttile
plot(record(:,1))
title("alpha="+n(1)+"")

nexttile
histogram(record(:,2))
title("GOR0 = "+x0(2)+"")
nexttile
plot(record(:,2))
title("alpha="+n(2)+"*current")

nexttile
histogram(record(:,3))
title("WC0 = "+x0(3)+"")
nexttile
plot(record(:,3))
title("alpha="+n(3)+"*current")

sgtitle(" reject rate = "+RejectRate/nIter+"") 
saveas(gcf,"Chain1.png")
close all;