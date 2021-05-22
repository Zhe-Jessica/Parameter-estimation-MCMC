%% scatter plot
%[PI, GOR, WC]
close all; clear all; clc;
load('Chain1.mat')
tiledlayout(1,3) 

nexttile
scatter(record(:,1),record(:,2))
xlabel('PI') 
ylabel('GOR') 

nexttile
scatter(record(:,1),record(:,3))
xlabel('PI') 
ylabel('WC') 

nexttile
scatter(record(:,2),record(:,3))
xlabel('GOR') 
ylabel('WC') 

% sgtitle("p_{PI}= "+n(1)+"          p_{GOR}= "+n(2)+"          p_{WC}= "+n(3)+"") 
