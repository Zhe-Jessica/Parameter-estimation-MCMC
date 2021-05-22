% Calculation of offline average and standard deviation matrixes  
close all; clc;
load('ym.mat')
%%  calculate average, standard deviation
% N = [440, 1440, 1620, 100];
% point = [0, 440, 1880, 3500, 3600];
% Ave = zeros(7,4);
% Cov = zeros(7,4);
%  for k = 2:5   
%      for j = 2:8
%          segment = ym(j, point(k-1)+1:point(k));
%          Ave(j-1,k) = mean(segment);
%          Cov(j-1,k) = std(segment);        
%      end 
%  end
%   Ave(:,1) = [];
%   Cov(:,1) = [];  
% %   save('cal.mat','Ave','Cov','N','point')
% save('offlinematrixes.mat','Ave','Cov','N','point')
% clear all; 

%% for one step input
 ym(:,1:1000) = []; 
 for j = 2:8
     segment = ym(j, :);
     Ave(j-1) = mean(segment);
     Cov(j-1) = std(segment);        
 end 
% Ave(:,1) = [];
% Cov(:,1) = [];  
%   save('cal.mat','Ave','Cov','N','point')
save('onematrixes.mat','Ave','Cov')
% clear all; 