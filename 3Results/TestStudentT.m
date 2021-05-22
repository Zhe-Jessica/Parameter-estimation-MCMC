% S = 10
% ye = 200
% N =1000;
% for i = 1:1000
% x = trnd(N -1)
% d = x*S/(N-1)^(-1/2)+ye
% plo(i) = d;
% end
load('Results.mat')
load('cal.mat')
for j =1
    for k = 1:4  
         for i = 1: N(k)
             for theta = 1:1000
            x = trnd(N(k) -1);
            d = x*Cov(j,k)/(N(k) -1)^(1/2) +Results(j,point(k)+i,theta);                   
            sample(j,point(k)+i,theta) = d;
             end
         end
    end
end
