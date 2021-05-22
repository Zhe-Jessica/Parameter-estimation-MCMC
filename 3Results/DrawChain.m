close all; clear all; clc;
for chain = 1:1
for i = 1:2   
%% result
load("test"+string(chain)+"_"+string(i)+".mat")
% load('test'+string(i)+'.mat')
% plot histogram of three parameters + x0, sigma 
figure
tiledlayout(3,2)
nexttile
histogram(record(:,1))
title("PI0 = "+x0(1)+"")
nexttile
plot(record(:,1))
title("sigma="+sigma(i)+"")

nexttile
histogram(record(:,2))
title("GOR0 = "+x0(2)+"")
nexttile
plot(record(:,2))
title("alpha="+n(i)+"*current")

nexttile
histogram(record(:,3))
title("WC0 = "+x0(3)+"")
nexttile
plot(record(:,3))
title("alpha="+n(i)+"*current")

sgtitle(" reject rate = "+RejectRate/nIter+"") 
end
end