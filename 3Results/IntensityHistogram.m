%% 2D-intensity histogram 
%[PI, GOR, WC]
close all; clear all; clc;
load('90000.mat');
tiledlayout(1,3) 

nexttile
X1 = [record(:,1),record(:,2)];
hist3(X1)
xlabel('PI')
ylabel('GOR')

nexttile
X2 = [record(:,1),record(:,3)];
hist3(X2)
xlabel('PI') 
ylabel('WC') 

nexttile
X3 = [record(:,2),record(:,3)];
hist3(X3)
xlabel('GOR') 
ylabel('WC') 

sgtitle("p_{PI}= "+n(1)+"          p_{GOR}= "+n(2)+"          p_{WC}= "+n(3)+"") 