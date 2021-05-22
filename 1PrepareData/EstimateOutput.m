function ye= EstimateOutput(PI, GOR, WC)
% load('meaurements.mat');
% ym = ym(:, 1:4600);
% save('ym.mat','ym')

 load('ym.mat');
 times = ym(1,:);
 inputs = ym(9,:);
 samples_number= length(times);
 ye = zeros(7,samples_number);
 
m0 =  [10728.107920;2878.73920;17222.344042]; 

 for k = 1:samples_number-1              % cut the begining of data
     tspan = [times(k) times(k+1)]; 
     u = inputs(k);                                      % apply input after Tdelay
     options = odeset('RelTol',1e-6,'AbsTol',1e-10);
     [t,m] = ode15s(@(t,m)  GLOWmodel( t , m , u, PI, GOR, WC), tspan, m0, options);
     m0 =  m(end,:)';     
     y = FindOtherStates(t, m, PI, GOR, WC);
     ye(:,k+1) = y(end,:)'; 
 end
%  figure
  ye(:,1:1000) = []; 
%   ym(:,1:1000) = []; 
  times(1:1000) = []; 
% plot(times,ye)
end