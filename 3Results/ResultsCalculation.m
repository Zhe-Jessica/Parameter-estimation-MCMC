% Calculate and save Results of the outputs of model with a chains of
% parameters
close all; clear all; clc;
load('90000.mat');
load('ym.mat');
times = ym(1,:);
inputs = ym(9,:);
samples_number= length(times);
ye = zeros(7,samples_number);

for iter = 1:1000
    iter
     m0 =  [10728.107920;2878.73920;17222.344042]; 
     theta =record(iter,:);
     PI = theta(1);
     GOR = theta(2);
     WC = theta(3);
     for k = 1:samples_number-1              % cut the begining of data
         tspan = [times(k) times(k+1)]; 
%          u = inputs(k);                                      % apply input after Tdelay
         u = 50; 
         options = odeset('RelTol',1e-6,'AbsTol',1e-10);
         [t,m] = ode15s(@(t,m)  GLOWmodel( t , m , u, PI, GOR, WC), tspan, m0, options);
         m0 =  m(end,:)';     
         y = FindOtherStates(t, m, PI, GOR, WC);
         ye(:,k+1) = y(end,:)'; 
     end
     ye(:,1:1000) = []; 
     Results(:,:,iter) = ye;
end
save('Results.mat','Results');
% times(1:1000) = []; 
% ym(:,1:1000) = []; 
%  figure

% plot(times,ye)
