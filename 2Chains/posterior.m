% need offline matrixes
function [P,ye]= posterior (PI, GOR, WC)
load('onematrixes.mat')
ye = EstimateOutput(PI, GOR, WC);
P = 0;
% count = 0;
%% multi step inpts
%     for j =1:7
%         addk = 0;
%         for k = 1:4  
%             inner = 0; 
%              for i = 1: N(k)
%                 inner = log(1+(Ave(j,k)-ye(j,point(k)+i))^2/Cov(j,k)^2)+ inner;
%     %             count = count +1;
%              end
%              inner = -N(k) * inner/2;
%              addk = addk+inner;
%         end
%         P = P+addk;
%     end
    
%% one step input
    for j =1:7
         inner = 0; 
         for i = 1: 3600
            inner = log(1+(Ave(j)-ye(j,i))^2/Cov(j)^2)+ inner;
         end
         inner = -3600 * inner/2;
         P = P+inner;
    end
end