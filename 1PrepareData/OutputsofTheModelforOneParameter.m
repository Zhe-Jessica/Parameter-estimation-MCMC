% Calculate the output of model with one theta and input
close all; clear all; clc;
%% Define the constants
load('data.mat');                                               % include inputs and time
 times = data(1,:);
 inputs = data(2,:);
 samples_number= length(times);
 ye = zeros(7,samples_number);
 ye(:,1) = data(3:9,1)';                                     % initial values of outputs
m0 =  [10728.107920;2878.73920;17222.344042]; 

% theta =[99999.9713, 0.1825, 0.1651];       %[PI, GOR, WC]
theta =[96000, 0.9, 0.94]
PI = theta(1);
GOR = theta(2);
WC = theta(3);
     
%% y_e                                                  estimated output with paramter
 for k = 1:samples_number-1;              % cut the begining of data
     tspan = [times(k) times(k+1)]; 
     u = inputs(k);                                      % apply input after Tdelay
     options = odeset('RelTol',1e-6,'AbsTol',1e-10);
     [t,m] = ode15s(@(t,m)  GLOWmodel( t , m , u, PI, GOR, WC), tspan, m0, options);
     m0 =  m(end,:)';     
     y = FindOtherStates(t, m, PI, GOR, WC);
     ye(:,k+1) = y(end,:)'; 
 end
 
 %% present result, compare outputs of model and measurements
% plot(times,ye,'-o')
% % hold on 
% plot(times,data(3:9,:)','-o')
figure
plot(times,ye(5:7,:),'--',times,data(7:9,:)')
figure
plot(times,ye(1:4,:),'--',times,data(3:6,:)')
% plot(times(500:8923),ye(1,500:8923),'--',times(500:8923),data(3,500:8923)')

