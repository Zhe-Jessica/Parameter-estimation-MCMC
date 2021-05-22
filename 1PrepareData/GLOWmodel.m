function dx = GLOWmodel( t , m , u, PI, GOR, WC)
%GLOWmodel: Summary of this function goes here
%   This is a model of one gas lifted oil well. X(dim=4*1) includes 4
%   differential variables (states) as belows:
%   X=[m_gp(t), m_ga(t), m_gt(t), m_lt(t)]'
%   DIM(theta): 1*4
%   theta = [WC PI GOR wn]    wn=well number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define the constants                                                  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Valve constants for gas lift choke valve (first element) and production
% % choke valve (the second element). 
valve_constants = [27.3, 273];
N6 = valve_constants(1,1); Nb6 = valve_constants(1,2);
Z_pc = 0.7076;
P_c = 200;
%%constat used in gas expansion formula
aY = 0.66;
% % universal gas constant
R = 8.314e-5;                                                   %% [m3.bar/K.mol]
% % fixed temperature assumed for gas everywhere in system
T = 280;                                                        %% [K]
% % molar mass of the lift gas
M = 20e-3;                                                      %% [Kg/mol]
% % gravitational acceleration02]
g = 9.81;                                                       %% [m/s2]

% % Water cut for five wells in order. DIM(WC):1*5
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
% L_pt = 13000;                                                   %% [m]
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
P_s=30; %% pressure of the common gathering manifold            %% [bar]
P_r = 150; %% pressure of reservoir       %% [bar]
% % minimum pressures in the gas distribution pipeline, in the annulus at
% % the point of injection, and in the tubing at the well head.
% % P_min = [P_cmin, P_ainjmin, P_whmin]. DIM(P_min):1*3
P_min = [10, 10, 10];                                           %% [bar]   %% ***(my assumption)***
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define initial states                                                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % extracting the states from X vector.
% %  mass of the gas in the distribution pipeline
%m_gp = X(1);                                                    %% [kg]
% %  mass of the gas in the annulus
m_ga = m(1);                                             %% [kg]
% %  mass of the gas in the tubing above the injection point
m_gt = m(2);                                            %% [kg]
% %  mass of the liquid (oil and water) in the tubing above the injection point
m_lt = m(3);                                           %% [kg]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define the control inputs                                             %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % supplied lift gas
w_gc = 8000;                                                   %% [Sm3/hr]
% % % in case we have step change in w_gc
% if t>25
%     w_gc = 6000;                                             %% [Sm3/hr]
% end
w_gc = w_gc*0.68; %% converting Sm3/hr to kg/hr                 %% [kg/hr]

% % valve opening of the gas lift choke valve
% u1 = 80;
% % % in case we have step change in u1
% if t>25
%     u1=30;
% end

% % valve opening of the production choke valve
u2 = 100;
% % % in case we have step change in u2
% if t>25
%     u2=55;
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define all the other variables                                              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Algebaic variables
% % Valve characteristic as a function of its opening for the lift gas
% % choke valve and the production choke valve. EQ(1) & EQ(2)
Cv1 = 0; Cv2 = 0; %% initialization
if (5<u) && (u<50)
Cv1 = 0.111*u-0.556;
elseif u>=50
Cv1 = 0.5*u-20;
end

if (5<u2) && (u2<50)
Cv2 = 0.111*u2-0.556;
elseif u2>=50
Cv2 = 0.5*u2-20;
end

% % Density of gas in the gas distribution pipeline: EQ(3)
%rho_gp = m_gp/(A_p*L_pt);                                       %% [Kg/m3] (3)
rho_gp = P_c*M/(Z_pc * R * T);

% % Density of gas in the annulus: EQ(4)
rho_ga = m_ga / (A_a*L_at);                                   %% [Kg/m3]

% % The average density of the mixture of liquid and gas in the tubing
% % above the injection point: EQ(5)
rho_m = (m_gt+m_lt) / (A_t*L_tt);                             %% [Kg/m3]

% % Z_pc = -2.572e-8*(P_c^3)+2.322e-5*(P_c^2)-0.005077*P_c+1; EQ(6)
% But, I assumed constant Z at 200 bar


% % Pressure of gas in the gas distribution pipeline upstream the lift gas
% % choke valve: This is scaler, not a vector. EQ(7)
%P_c = Z_pc*m_gp*R*T/(M*A_p*L_pt);                               %% [bar]


% % Z_pa1 = -2.572e-8*(P_a1^3)+2.322e-5*(P_a1^2)-0.005077*P_a1+1; EQ(8)
% But, I assumed constant Z at 170 bar.
Z_pa = 0.6816;

% % Pressure of gas in the annulus downstream the lift gas choke valve: EQ(9)
P_a = (R*T/M) * (Z_pa*m_ga) / (A_a*L_at);                    %% [bar]

% % Pressure of gas in the annulus upstream the gas injection valve. EQ(10) 
P_ainj = P_a + 1e-5*g*(rho_ga*L_av);                           %% [bar]

% % Gas expansion factor in the gas lift choke valve. EQ(11) 
Y1 = 1 - aY*(P_c-P_a) / max(P_c,P_min(1,1));

% % Density of liquid (oil and water) in each well based on its water cut. EQ(12) 
rho_l = Rho(1,2)*WC + Rho(1,1)*(1-WC);                  %% [Kg/m3]

% % The volume of gas present in the tubing above the gas injection point. EQ(13)
V_G = (A_t*L_tt) - (m_lt/rho_l);                             %% [m3]

% % Z_PG(i) = -2.572e-8*(PG(i)^3)+2.322e-5*(PG(i)^2)-0.005077*PG(i)+1;
% % where PG(i)=(P_wh(i)+P_tinj(i))/2. But, I assumed constant Z at 150 bar
% % EQ(14)
Z_PG = 0.6741;

% % dP_ft1 = 0.5*f_Dt1*L_tt1*rho_m1*(v_t1^2)*D_h1; Pressure loss due to
% % friction from the injection point to well head. But, I assumed no
% % friction. EQ(15)
dP_ft = 0;                                        %% [bar]

% % Pressure in the tubing downstream the gas injection valve. EQ(16)
P_tinj = (R*T/M)*((Z_PG*m_gt)/V_G) + ...
         0.5*g*1e-5*(rho_m*L_tv) + 0.5*dP_ft;                  %% [bar]

% % Pressure in the tubing upstream the production choke valve. EQ(17)
P_wh = (R*T/M)*((Z_PG*m_gt)/V_G) - ...
         0.5*g*1e-5*(rho_m*L_tv) - 0.5*dP_ft;                  %% [bar]

% % dP_fr1 = 0.5*f_Dr1*L_rt1*rho_r*(v_r1^2)*D_h1; Pressure loss due to
% % friction from the bottom hole to the injection point. But, I assumed no
% % friction. DIM(dP_fr):1*5. EQ(18)
dP_fr = 0;                                        %% [bar]

% % The bottom hole pressure or well flow pressure. EQ(19)
P_wf = P_tinj + (rho_r*g*1e-5*L_rv) + dP_fr;                    %% [bar]

% % Gas expansion factor in the gas injection valve. EQ(20)
Y2 = 1 - aY*(P_ainj-P_tinj)/max(P_ainj,P_min(1,2));

% % Gas expansion factor in the production choke valve. EQ(21)
Y3 = 1 - aY*(P_wh-P_s)/max(P_wh,P_min(1,3));

% % MASS FLOW RATES
% % Mass flow rate of the gas through the gas lift choke valve. EQ(22)
w_ga = N6*Cv1*Y1*sqrt(rho_gp*max(P_c-P_a,0)); %% [Kg/hr]
    
% % Mass flow rate of the gas injected into the tubing from the annulus. EQ(23)
w_ginj = K*Y2*sqrt(rho_ga*max(P_ainj-P_tinj,0)); %% [Kg/hr]

% % The mass flow rate of the liquid from the reservoir. EQ(24)
w_lr = PI*max(P_r-P_wf,0); %% [Kg/hr]
    
% % The mass flow rate of gas from the reservoir. EQ(25)
w_gr = GOR*w_lr; %% [Kg/hr]

% % The mass flow rate of the mixture of gas and liquid from the production
% % choke valve. EQ(26)
w_gop = Nb6*Cv2*Y3*sqrt(rho_m*max(P_wh-P_s,0)); %% [Kg/hr]

% % Mass flow rate of gas through the production choke valve. EQ(27)
w_gp = m_gt*w_gop/(m_gt+m_lt); %% [Kg/hr]

% % Mass flow rate of liquid through the production choke valve. EQ(28)
w_lp = m_lt*w_gop/(m_gt+m_lt); %% [Kg/hr]

% % Oil compartment mass flow rate from liquid product considering water cut. EQ(29)
w_op = (1-WC)*w_lp*Rho(1,1)/rho_l; %% [Kg/hr]
    
% % Water compartment mass flow rate from liquid product considering water cut. EQ(30)
w_wp = w_lp - w_op; %% [Kg/hr]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% % Define the derivatives                                              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Mass balance in gas distribution manifold:
% % This is just a single scaler, not a vector. EQ(31)
%dm_gpdt = w_gc - w_ga; %% [Kg/hr]

% % Mass balance in annulus. EQ(32)
dm_gadt = w_ga - w_ginj; %% [Kg/hr]

% % Mass balance for the gas in tubing and above the injection point. EQ(33)
dm_gtdt = w_ginj + w_gr - w_gp; %% [Kg/hr]

% % Mass balance for the liquid in tubing and above the injection point. EQ(34)
dm_ltdt = w_lr - w_lp; %% [Kg/hr]

dx = [dm_gadt; dm_gtdt; dm_ltdt];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%