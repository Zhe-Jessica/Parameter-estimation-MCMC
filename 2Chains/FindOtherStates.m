function y = FindOtherStates(t, x, PI, GOR, WC)
% Summary of this function goes here
%   This is a model of five gas lifted oil wells. X(dim=16*1) includes 16
%   differential variables as belows:
%   X=[m_gp(t), m_ga(t), m_gt(t), m_lt(t)]'
%   DIM(theta): 1*3
%   theta = [WC PI GOR wn]    wn=well number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define the constants                                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Valve constants for gas lift choke valve (first element) and production
% % choke valve (the second element). 
valve_constants = [27.3, 273];
N6 = valve_constants(1,1); Nb6 = valve_constants(1,2);

%%constat used in gas expansion formula
aY = 0.66;
% % universal gas constant
R = 8.314e-5;                                                   %% [m3.bar/K.mol]
% % fixed temperature assumed for gas everywhere in system
T = 280;                                                            %% [K]
% % molar mass of the lift gas
M = 20e-3;                                                       %% [Kg/mol]
% % gravitational acceleration02]
g = 9.81;                                                           %% [m/s2]

%3% Water cut for five wells in order. DIM(WC):1*5
% % WC = [WC1, WC2, WC3, WC4, WC5]
% WC = [0.25, 0.04, 0.15, 0.16, 0.15];
% WC = theta(1);

% % number of well
% nw = theta(4);
nw = 1;

% % gas injection valve constants in order of well numbers DIM(K):1*5
Kws = [68.43, 67.82, 67.82, 69.26, 66.22];                        %% [(Kg.m3/bar)^.5/hr]
K = Kws(nw);

% % Productivity Index in order of well numbers. DIM(PI):1*5
% PI = [2.51e4, 1.63e4, 1.62e4, 1.75e4, 2.232e4];                 %% [kg/hr.bar]
% PI = theta(2);

% % Gas to Oil Ratio. DIM(GOR):1*5
% GOR = [0.10, 0.15, 0.08, 0.10, 0.06];                                      %% ***(my assumption)***
% GOR = theta(3);

% % Densities of crude oil, water, and fluid mixture(oil, water, and gas)
% % Rho = [rho_o, rho_w, rho_r]
Rho = [800, 1000, 700];                                         %% [Kg/m3]
rho_r = Rho(1,3);                                               %% [Kg/m3]

% % PIPE DIMENSIONS
% % DIAMETERS
% % Lift gas distribution pipeline
IDp = 3;                                                        %% [in]    ***(my assumption)***
% %  Tubing inner and outer diameters
IDtws = [6.18, 6.18, 6.18, 6.18, 6.18];                           %% [in]
IDt = IDtws(nw);
ODtws = [7.64, 7.64, 7.64, 7.64, 7.64];                           %% [in]
ODt = ODtws(nw);
% % Annulus diameters
IDaws = [9.63, 9.63, 9.63, 9.63, 9.63];                           %% [in]
IDa = IDaws(nw);

% % Length
% % Lift gas distribution pipeline diameter and length
L_pt = 13000;                                                   %% [m]
% % total length of tubings above the injection point
L_ttws = [2758, 2559, 2677, 2382, 2454];                          %% [m]
L_tt = L_ttws(nw);
% % vertical length of tubings above the injection point
L_tvws = [2271, 2344, 1863, 1793, 1789];                          %% [m]
L_tv = L_tvws(nw);
% % vertical length of tubings below the injection point
L_rvws = [114, 67, 61, 97, 146];                                  %% [m]
L_rv = L_rvws(nw);
% % total length of annulus
L_at  = L_tt;                                                   %% [m]
% % vertical length of annulus
L_av = L_tv;                                                    %% [m]

% % Cross section areas
% % gas distribution pipeline
A_p = pi*((IDp*0.0254)^2)/4;                                    %% [m2]
% % Tubing
A_t = pi*((IDt*0.0254)^2)/4;                                   %% [m2]
% % Annulus
A_a = (pi/4)*((IDa*0.0254)^2-(ODt*0.0254)^2);                 %% [m2]

% % pressures all in [bar]
P_r = 150; %% pressure of reservoir       %% [bar]
% % minimum pressures in the gas distribution pipeline, in the annulus at
% % the point of injection, and in the tubing at the well head.
% % P_min = [P_cmin, P_ainjmin, P_whmin]. DIM(P_min):1*3
P_min = [10, 10, 10];                                           %% [bar]   %% ***(my assumption)***

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define initial states                                                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = size(t,1); %% time length
% m_gp = X(:,1);
m_ga=x(:,1);
m_gt=x(:,2);
m_lt=x(:,3);
% % preallocating the variable to increase the speed
Cv1 = zeros(N,1);
Cv2 = zeros(N,1);
rho_gp = zeros(N,1);
rho_ga = zeros(N,1);
rho_m = zeros(N,1);
rho_l = zeros(N,1);
Z_pc = zeros(N,1); 
Z_pa = zeros(N,1);
P_a = zeros(N,1);
P_ainj = zeros(N,1);
Y1 = zeros(N,1);
Y2 = zeros(N,1);
Y3 = zeros(N,1);
V_G = zeros(N,1);
Z_PG = zeros(N,1);
dP_ft = zeros(N,1);
P_s= zeros(N,1); %% pressure of the common gathering manifold            %% [bar]
P_c= zeros(N,1); 
P_tinj = zeros(N,1);
P_wh = zeros(N,1);
dP_fr = zeros(N,1);
P_wf = zeros(N,1);
w_ga = zeros(N,1);
w_ginj = zeros(N,1);
w_lr = zeros(N,1);
w_gr = zeros(N,1);
w_gop = zeros(N,1);
w_gp = zeros(N,1);
w_lp = zeros(N,1);
w_op = zeros(N,1);
w_wp = zeros(N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define the control inputs for whole the simulation time               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % supplied lift gas
% w_gc = 8000;                                                    %% [Sm3/hr]
% w_gc = w_gc*0.68*ones(N,1); %% converting Sm3/hr to kg/hr        %% [kg/hr]
% % % valve opening of the gas lift choke valves
u = 80;
u = ones(N,1)*u;
% % % valve opening of the production choke valves
u2 = 100;
u2 = ones(N,1)*u2;

%% Define all the variables
% % algebaic variables
for i=1:1:N
    
    % % % in case we have step change in control inputs
%     if t(i,1)>25
%         w_gc(i,1)=6000*0.68;
%     end
%     if t(i,1)>25
%         u1(i,1)=30;
%     end

    % % Valve characteristic as a function of its opening for the lift gas
    % % choke valves and the production choke valves. EQ(1) & EQ(2)

    if (5<u(i)) && (u(i)<50)
    Cv1(i) = 0.111*u(i)-0.556;
    elseif u(i)>=50
    Cv1(i) = 0.5*u(i)-20;
    end

    if (5<u2(i)) && (u2(i)<50)
    Cv2(i) = 0.111*u2(i)-0.556;
    elseif u2(i)>=50
    Cv2(i) = 0.5*u2(i)-20;
    end


    
    % % Density of gas in the annulus: DIM(rho_ga):N*1. EQ(4)
    rho_ga(i) = m_ga(i) / (A_a*L_at);

    % % The average density of the mixture of liquid and gas in the tubing
    % % above the injection point: DIM(rho_m):N*1 EQ(5)
    rho_m(i) = (m_gt(i)+m_lt(i)) / (A_t*L_tt);
    
    % % Z_pc = -2.572e-8*(P_c^3)+2.322e-5*(P_c^2)-0.005077*P_c+1; EQ(6)
    % But, I assumed constant Z at 200 bar
    Z_pc(i) = 0.7076;
    
    % % Pressure of gas in the gas distribution pipeline upstream the lift gas
    % % choke valve: DIM(P_c):N*1. EQ(7)
%     P_c(i) = Z_pc(i)*m_gp(i)*R*T/(M*A_p*L_pt);            %% [bar]
    P_c(i) = 200;
    P_s(i) =30;
    
        % % Density of gas in the gas distribution pipeline: DIM(rho_gp):N*1. EQ(3)
%     rho_gp(i) = m_gp(i)/(A_p*L_pt);
    rho_gp(i) = P_c(i)*M/(Z_pc(i) * R * T);
    
    % % Z_pa1 = -2.572e-8*(P_a1^3)+2.322e-5*(P_a1^2)-0.005077*P_a1+1; EQ(8)
    % But, I assumed constant Z at 170 bar. DIM(Z_pa):N*1
    Z_pa(i) = 0.6816;
    
    % % Pressure of gas in the annulus downstream the lift gas choke valve:
    % % DIM(P_a):N*1. EQ(9)
    P_a(i) = (R*T/M) * (Z_pa(i)*m_ga(i)) / (A_a*L_at); %% [bar]
    
    % % Pressure of gas in the annulus upstream the gas injection valve:
    % % DIM(P_ainj):N*1. EQ(10) 
    P_ainj(i) = P_a(i) + 1e-5*g*(rho_ga(i)*L_av);        %% [bar]

    % % Gas expansion factor in the gas lift choke valve:
    % % DIM(Y1):N*1. EQ(11) 
    Y1(i) = 1 - aY*(P_c(i)-P_a(i)) / max(P_c(i),P_min(1,1));
    
    % % Density of liquid (oil and water) in each well based on its water cut:
    % % DIM(rho_l):N*1. EQ(12)
    rho_l(i) = Rho(1,2)*WC + Rho(1,1)*(1-WC);         %% [Kg/m3]
    
    % % The volume of gas present in the tubing above the gas injection point:
    % % DIM(rho_l):N*1. EQ(13)
    V_G(i) = (A_t*L_tt) - (m_lt(i)/rho_l(i)); %% [m3]
    
    % % Z_PG(i) = -2.572e-8*(PG(i)^3)+2.322e-5*(PG(i)^2)-0.005077*PG(i)+1;
    % % where PG(i)=(P_wh(i)+P_tinj(i))/2. But, I assumed constant Z at 150 bar
    % % DIM(Z_PG):N*1. EQ(14)
    Z_PG(i) = 0.6741;

    % % dP_ft1 = 0.5*f_Dt1*L_tt1*rho_m1*(v_t1^2)*D_h1; Pressure loss due to
    % % friction from the injection point to well head. But, I assumed no
    % % friction: DIM(dP_ft):N*1. EQ(15)
    dP_ft(i) = 0;                               %% [bar]
    
    % % Pressure in the tubing downstream the gas injection valve:
    % % DIM(P_tinj):N*1. EQ(16)
    P_tinj(i) = (R*T/M)*((Z_PG(i)*m_gt(i))/V_G(i)) + ...
                  0.5*g*1e-5*(rho_m(i)*L_tv) + 0.5*dP_ft(i); %% [bar]
              
    % % Pressure in the tubing upstream the production choke valve:
    % % DIM(P_wh):N*1. EQ(17)
    P_wh(i) = (R*T/M)*((Z_PG(i)*m_gt(i))/V_G(i)) - ...
                  0.5*g*1e-5*(rho_m(i)*L_tv) - 0.5*dP_ft(i); %% [bar]

    % % dP_fr1 = 0.5*f_Dr1*L_rt1*rho_r*(v_r1^2)*D_h1; Pressure loss due to
    % % friction from the bottom hole to the injection point. But, I assumed no
    % % friction. DIM(dP_fr):N*1. EQ(18)
    dP_fr(i) = 0; %% [bar]
    
    % % The bottom hole pressure or well flow pressure. EQ(19)
    P_wf(i) = P_tinj(i) + (rho_r*g*L_rv*1e-5) + dP_fr(i); %% [bar]
    
    % % Gas expansion factor in the gas injection valve: DIM(Y2):N*1. EQ(20)
    Y2(i) = 1 - aY*(P_ainj(i)-P_tinj(i))/max(P_ainj(i),P_min(1,2));

    % % Gas expansion factor in the production choke valve: DIM(Y3):N*1. EQ(21)
    Y3(i) = 1 - aY*(P_wh(i)-P_s(i))/max(P_wh(i),P_min(1,3));
    
    % % MASS FLOW RATES
    % % Mass flow rate of the gas through the gas lift choke valve:
    % % DIM(w_ga):N*1. EQ(22)
    w_ga(i) = N6*Cv1(i)*Y1(i)*sqrt(rho_gp(i)*max(P_c(i)-P_a(i),0)); %% [Kg/hr]
    
    % % Mass flow rate of the gas injected into the tubing from the annulus:
    % % DIM(w_ginj):N*1. EQ(23)
    w_ginj(i) = K*Y2(i)*sqrt(rho_ga(i)*max(P_ainj(i)-P_tinj(i),0)); %% [Kg/hr]
        
    % % The mass flow rate of the liquid from the reservoir:
    % % DIM(w_lr):N*1. EQ(24)
    w_lr(i) = PI*max(P_r-P_wf(i),0); %% [Kg/hr]
        
    % % The mass flow rate of gas from the reservoir:
    % % DIM(w_gr):N*1. EQ(25)
    w_gr(i) = GOR*w_lr(i); %% [Kg/hr]
        
    % % The mass flow rate of the mixture of gas and liquid from the
    %  % production choke valve: DIM(w_gop):N*1. EQ(26)
    w_gop(i) = Nb6*Cv2(i)*Y3(i)*sqrt(rho_m(i)*max(P_wh(i)-P_s(i),0)); %% [Kg/hr]
        
    % % Mass flow rate of gas through the production choke valve:
    % % DIM(w_gp):N*1. EQ(27)
    w_gp(i) = m_gt(i)*w_gop(i)/(m_gt(i)+m_lt(i)); %% [Kg/hr]
        
    % % Mass flow rate of liquid through the production choke valve:
    % % DIM(w_lp):N*1. EQ(28)
    w_lp(i) = m_lt(i)*w_gop(i)/(m_gt(i)+m_lt(i)); %% [Kg/hr]

    % % Oil compartment mass flow rate from liquid product considering water cut:
    % % DIM(w_op):N*1. EQ(29)
    w_op(i) = (1-WC)*w_lp(i)*Rho(1,1)/rho_l(i); %% [Kg/hr]
    
    % % Water compartment mass flow rate from liquid product considering water cut:
    % % DIM(w_wp):N*1. EQ(30)
    w_wp(i) = w_lp(i) - w_op(i); %% [Kg/hr]    
       
end
%% Define outputs1
% control_inputs = [w_gc, u1, u2];
%y = [w_ga, w_gp, w_op, w_wp, P_wf, P_wh, P_s, P_c, P_a];
y = [w_ga, w_gp, w_op, w_wp, P_wf, P_wh, P_a];
% pressures = [ , P_ainj, P_tinj];
% flows = [w_ginj, w_lr, w_gr,  w_lp, w_gop, ];
% den = [rho_m, rho_l];
end