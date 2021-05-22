% clear all; 
close all; clc;clear;
% offline matrixes
load('Results.mat')

% load("Chain"+chain+".mat")
load('onematrixes.mat')
%% posterior predictive using student t
for j =1:7
         for i = 1: 3600
             for theta = 1:1000
            x = trnd(3600 -1);
            d = x*Cov(j)/(3600 -1)^(1/2)+Results(j,i,theta);                   
            sample(j,i,theta) = d;
             end
         end
end

% FirstOutput = squeeze(sample(3,:,:));         %for test
    
%% quantile calcualtion and record
% median( A ) 
% quantile(x,4)
% quantile(x,0.30)

%  load('sample.mat')
Range =zeros(7,3);
for j =1:7
     for i = 1: 3600
       p_d (j,i,1:3)=quantile(sample(j,i,:),3)';
%        p_d (i,1:3)=quantile(FirstOutput(i,:),3)';
     end
    RangeLow = median(p_d(j,:,1)) ;
    RangeMed = median(p_d(j,:,2)) ;
    RangeHigh = median(p_d(j,:,3)) ;
    Range(j,1) = RangeLow;
    Range(j,2) = RangeMed;
    Range(j,3) = RangeHigh;
end
save('Pos3','sample', 'Range')

%% plot 
load('ym.mat', 'ym')
% ym(:,1:1000) =[];
figure
tiledlayout(1,3)%(2,2)
VectorA = zeros(1,4600);
for i = 5:7%1:4 
    RangeLow = ones(size(VectorA))* Range(i,1);
    RangeMed = ones(size(VectorA))*Range(i,2);
    RangeHigh = ones(size(VectorA))* Range(i,3);

    nexttile
     p(1).LineWidth = 5;
     p(2).LineWidth = 10;
    plot(ym(1,:), ym(i+1,:),'Color',[0 0 0]+0.05*15)
    xlim([0 25.6])
    hold on    
    plot(ym(1,:), RangeLow,'LineWidth',3)    
    hold on
    plot(ym(1,:), RangeMed,'LineWidth',3)  
    hold on
    plot(ym(1,:), RangeHigh,'LineWidth',3)    
    hold off
    xlabel('Time[hr]') 
    ylabel('Pressure[bar]') 
%     ylabel('Flow rate[kh/hr]') 
    title("Output"+i+"")
%     saveas(gcf,"Chain"+chain+".png")
%     close all
end