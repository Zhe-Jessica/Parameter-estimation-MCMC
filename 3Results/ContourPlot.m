close all; clc; clear;
load('Chain1.mat')
rng(5)
colormap jet
tiledlayout(3,2) 

%%
nexttile
scatter(record(:,1),record(:,2))
xlabel('PI') 
ylabel('GOR') 
title('PI vs. GOR')

nexttile
X=[record(:,1),record(:,2)];
gridx1 = linspace(min(X(:,1)),max(X(:,1)),500);
gridx2 = linspace(min(X(:,2)),max(X(:,2)),500);
% gridx1 = linspace(9.6e+04,max(X(:,1)),100);
% gridx2 = linspace(0.924,0.932,100);
[x1,x2] = meshgrid(gridx1, gridx2);
x1 = x1(:);
x2 = x2(:);
xi = [x1 x2];
ksdensity(X,xi,'PlotFcn','contour');
xlabel('PI') 
ylabel('GOR') 
colorbar()

%%
nexttile
scatter(record(:,1),record(:,3))
xlabel('PI') 
ylabel('WC') 
title('PI vs. WC')

nexttile
X=[record(:,1),record(:,3)];
gridx1 = linspace(min(X(:,1)),max(X(:,1)),500);
gridx2 = linspace(min(X(:,2)),max(X(:,2)),500);
% gridx1 = linspace(7.3e+04,10e+04,100);
% gridx2 = linspace(0.798,0.802,100);
[x1,x2] = meshgrid(gridx1, gridx2);
x1 = x1(:);
x2 = x2(:);
xi = [x1 x2];
ksdensity(X,xi,'PlotFcn','contour');
xlabel('PI') 
ylabel('WC') 
colorbar()

%%
nexttile
scatter(record(:,2),record(:,3))
xlabel('GOR') 
ylabel('WC') 
title('GOR vs. WC')

nexttile
figure
X=[record(:,2),record(:,3)];
gridx1 = linspace(min(X(:,1)),max(X(:,1)),500);
gridx2 = linspace(min(X(:,2)),max(X(:,2)),500);
[x1,x2] = meshgrid(gridx1, gridx2);
x1 = x1(:);
x2 = x2(:);
xi = [x1 x2];
ksdensity(X,xi,'PlotFcn','contour');
xlabel('GOR') 
ylabel('WC') 
colorbar()