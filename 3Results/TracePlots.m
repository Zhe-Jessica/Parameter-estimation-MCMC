close all; clc; clear
figure
tiledlayout(3,2)

%% PI
nexttile
for i = 1:5
    load("Chain"+i+".mat")
    plot(record(:,1))
    hold on 
     clear 
end
xlabel('Iterations') 
ylabel('Parameter values') 
title("PI")
legend('Location','southeast')
legend('Chain1','Chain2','Chain3','Chain4','Chain5')

clear
load('Chain1.mat')
nexttile
histogram(record(:,1))

% GOR
nexttile
for i = 1:5
    load("Chain"+i+".mat")
    plot(record(:,2))
    hold on 
     clear 
end
xlabel('Iterations') 
legend('Location','southeast')
legend('Chain1','Chain2','Chain3','Chain4','Chain5')
title("GOR")
clear
load('Chain1.mat')
nexttile
histogram(record(:,2))

% WC
nexttile
for i = 1:5
    load("Chain"+i+".mat")
    plot(record(:,3))
    hold on 
     clear 
end
xlabel('Iterations') 
legend('Location','southeast')
legend('Chain1','Chain2','Chain3','Chain4','Chain5')
title("WC")

clear
load('Chain1.mat')
nexttile
histogram(record(:,3))


% sgtitle("Trace plots of thress parameters") 