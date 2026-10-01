function [sys,par] = ptfe(par)

% industrial batch polymerization reactor PTFE
%
% Written by: Dr. Chenchen Zhou, July. 2024 ZJU

%u=[FA,FB,Fj,Tjin] manipulated variable
%y=[C_R,C_A,C_B,C_C,C_E,Tr,NgA,P,Tj,RA,Mn,Mw,PDI] Measured variable

import casadi.*

Nlambda0=MX.sym('Nlambda0',1); %  moment of active polymer
Nlambda1=MX.sym('Nlambda1',1); %  moment of active polymer
Nlambda2=MX.sym('Nlambda2',1); %  moment of active polymer
Nmu0=MX.sym('Nmu0',1); %  moment of dead polymer
Nmu1=MX.sym('Nmu1',1); %  moment of dead polymer
Nmu2=MX.sym('Nmu2',1); %  moment of dead polymer
N_R=MX.sym('N_R',1); % mol/L Concentration of free radical
N_A=MX.sym('N_A',1); % mol/L Concentration of A
N_B=MX.sym('N_B',1); % mol/L Concentration of B
N_C=MX.sym('N_C',1); % mol/L Concentration of C
N_E=MX.sym('N_E',1); % mol/L Concentration of E
Tr=MX.sym('Tr',1); % K Temperature of reactor
NgA=MX.sym('NgA',1); % mol Amount of A in gas phase of reactor
P=MX.sym('P',1); % Pa Reactor gas pressure
Tjacket=MX.sym('Tjacket',1); % K Temperature of reactor jacket water
MAin=MX.sym('MAin',1); %kg Accumulated A feed weight
MBin=MX.sym('MBin',1); %kg Accumulated B feed weight
MA_l=MX.sym('MA_l',1); %kg Accumulated mass of A dissolved into the liquid phase
Vl=MX.sym('Vl',1); %kg Accumulated mass of A dissolved into the liquid phase
% Tab=MX.sym('Tab',1); %K

% Input struct (optimization variables):
FA=MX.sym('FA',1); % kg/s Feed flow of A
FB=MX.sym('FB',1); % kg/s Feed flow of B
% Fjacket=MX.sym('Fjacket',1); % kg/s Flow rate of jacket water
% Fjacket_h=MX.sym('Fjacket_h',1); % kg/s Flow rate of jacket water
% Tjin=MX.sym('Tjin',1); % K Temperature of reactor jacket water inlet
alpha=MX.sym('alpha',1); % [-1,0] cold water (0,1] hot water


%  Uncertain parameters:
deltaH=MX.sym('deltaH',1); %172700; %J/mol Standard reaction heat
UA=MX.sym('UA',1); %16870; % J/K Heat-transfer coefficient * Heat-transfer area
Thot=MX.sym('Thot',1);%363; %K Maximum temperature of jacket out cooling water
Tcold=MX.sym('Tcold',1);%283; %K Minimum temperature of jacket out cooling water

%Certain parameters
Fmax_jh=1; %kg/s Maximum feed rate of hot jecket water
Fmin_jh=0; %kg/s Minimum feed rate of hot jecket water
Fmax_jc=10;30; %kg/s Maximum feed rate of cold jecket water
Fmin_jc=0; %kg/s Minimum feed rate of cold jecket water
Mmax_A=3250; %kg Maximum feed mass of A
Mlim_B=0.0020; %kg Feeding APS mass limit for changing the B feed rate
% Thot=363; %K Maximum temperature of jacket out cooling water
% Tcold=283; %K Minimum temperature of jacket out cooling water

% Modified Split-Range Logic (Parameterized)
alpha_switch = 0.0;

% 1. Heating Flow (Fjacket_h)
% Active for alpha > alpha_switch.
% Range [alpha_switch, 1] maps to [0, 100%]
u_h = (alpha - alpha_switch) / (1 - alpha_switch);
Fjacket_h = (Fmax_jh - Fmin_jh) * if_else(alpha > alpha_switch, u_h, 0);

% 2. Cooling Flow (Fjacket)
% Active for alpha < alpha_switch.
% Range [-1, alpha_switch] maps to [100%, 0%]
u_c = (alpha_switch - alpha) / (alpha_switch + 1);

% Clamp to [0, 1] and apply
F_frac_cold = if_else(alpha < alpha_switch, min(1, max(0, u_c)), 0);

Fjacket = (Fmax_jc - Fmin_jc) * F_frac_cold;

R=8.314; %J/(mol*K) Ideal gas constant
Vl0=4000;  %L Reactor liquid volumn
V=6000;  %L Reactor gas volumn
C_A0=0.22857145;
P0=1.6e5;
Tr0=340;
NgA0=P0*(V-Vl0)/(R*Tr0)/1000;

FA_max=0.694; %kg/s Maximum feed rate of A
FA_min=0; %kg/s Minimum feed rate of A
FB_max=0.0000004; %kg/s Maximum feed rate of B
FB_min=0; %kg/s Minimum feed rate of B

% Allow external override of parameters for Mismatch Simulation
if isfield(par, 'k10'), k10=par.k10; else, k10=1.13e16*1e1; end
if isfield(par, 'k20'), k20=par.k20; else, k20=3.62e13*1e2; end
if isfield(par, 'k70'), k70=par.k70; else, k70=3.38e4*1e5; end

k30=5.49e6*1e1;1.249e7;
k40=3.32e6*3e1;
k50=1051*5e2*6;
k60=1051*1e3*0;
% k70 handled above
E10=1.348e5; %J/mol Activation energy
E20=119715;
E30=17413.76;
E40=53020;
E50=20000;
E60=20000;
E70=13604.5;

k=0.0562; %kg*m3/(s*mol)  Gas liquid mass transfer rate constant of A (NgA*1e3)
H=1/700; % mol/(Pa*m3)  Henry constant HP=c
% Vl0=2000;  %L Reactor liquid volumn
% V=6000;  %L Reactor gas volumn
Mg_A=100; %g/mol Relative molecular mass of A
Mg_B=228; %g/mol Relative molecular mass of B
Cp_A=804; %J/(kg*K) Specific heat capacity of A
TA=333; %K Tempereture of A feed
% UA=16870; % J/K Heat-transfer coefficient * Heat-transfer area
% deltaH=172700; %J/mol Standard reaction heat
a=1; %coefficient of heat loss to the ambience
b=2; %coefficient of heat loss to the ambience
Tamb=293; %K Temperature of ambience
Qstir=0.5; %W Stirring heat
M=4000; %kg Mass of reactor liquid phase
Mj=2000; %kg Mass of jacket water
Cp=4200; %J/(kg*K) Specific heat capacity of reactor liquid phase
Cp_p = 1300; %J/(kg*K) 
% Mmax_A=3250; %kg Maximum feed mass of A

Trsp=351; %K Reactor temperature set point
Psp=1500000; %Pa Reactor gas pressure set point
% Pact_Tc=1000000; %Pa Pressure standard for starting the temperature control

% algebraic equations

lambda0=Nlambda0/Vl; %  moment of active polymer
lambda1=Nlambda1/Vl; %  moment of active polymer
lambda2=Nlambda2/Vl; %  moment of active polymer
mu0=Nmu0/Vl; %  moment of dead polymer
mu1=Nmu1/Vl; %  moment of dead polymer
mu2=Nmu2/Vl; %  moment of dead polymer
C_R=N_R/Vl; % mol/L Concentration of free radical
C_A=N_A/Vl; % mol/L Concentration of A
C_B=N_B/Vl; % mol/L Concentration of B
C_C=N_C/Vl; % mol/L Concentration of C
C_E=N_E/Vl; % mol/L Concentration of E




k1 = k10*exp(-E10/(R*Tr));
k2 = k20*exp(-E20/(R*Tr));
k3 = k30*exp(-E30/(R*Tr));
k4 = k40*exp(-E40/(R*Tr));
k5 = k50*exp(-E50/(R*Tr));
k6 = k60*exp(-E60/(R*Tr));
k7 = k70*exp(-E70/(R*Tr));


Fgl = k*(P*H-1000*C_A);
M_t  = M; %+ NgA0*(Mg_A/1000) + MBin;
% Vl=Vl0+MA_l*0.6667;
Vg=V-Vl;
f = 1;

X=(C_A0*Vl0+MA_l*(1000/Mg_A)-C_A*Vl)/(C_A0*Vl+MA_l*(1000/Mg_A));

% For 3 effects
Mwcr=6e4;
Vfcr=0.129;
Vfcr2 = exp(-0.7-1000/Tr);
Vfcr3 = 0.039; %0.069
TgM=167; % k
Tgp=383; % K
% K3=exp(1.4+7460/Tr);
% m=0.5;
% n=1.75;
pA=1.11;
pB=1;
pC=0.42;
epsilon=0.2170;
epsilon=(0.966471-1.164e-3*Tr)/(1.19504-3.3e-4*Tr);
ap = 4.8e-4;
aM = 1e-3;


% free-volume
phi_p=(1+epsilon)*X/(1+epsilon*X);
phi_M=1-phi_p;
Vf = 0.025+ap*(Tr-Tgp)*phi_p+aM*(Tr-TgM)*phi_M;
K3 = exp(-0.4+4460/Tr);
K = (Mg_A*(Nmu2+Nlambda2)/(Nmu1+Nlambda1))^0.5*exp(1.11/Vf);
% Glass effect
% k3=k3*exp(-pB*0.5*((1/Vf-1/Vfcr2)+sqrt((1/Vf-1/Vfcr2)^2+1e-5)));
% k3=k3*exp(-pB*max(1/Vf-1/Vfcr2,0));

% Cage effect
% f = f*exp(-pC*0.5*((1/Vf-1/Vfcr3)+sqrt((1/Vf-1/Vfcr3)^2+1e-5)));
% f = f0*exp(-pC*max(1/Vf-1/Vfcr3,0));
% Gel effect
% K = Mw^m*exp(pA/Vf)-K3;>0
if par.isFault==1
    % k5=k5*exp(-pA*0.5*((1/Vf-1/Vfcr)+sqrt((1/Vf-1/Vfcr)^2+1e-5)));
    % k6=k6*exp(-pA*0.5*((1/Vf-1/Vfcr)+sqrt((1/Vf-1/Vfcr)^2+1e-5)));
    % if K>=K3 
    %     k7=k7*exp(-pA*0.5*((1/Vf-1/Vfcr)+sqrt((1/Vf-1/Vfcr)^2+1e-5)));
    % end
    Mwcr_true = Mg_A * (Nmu2 + Nlambda2) / (Nmu1 + Nlambda1);
    Mwcr_r = MX.sym('Mwcr_r');
    tau = 1e-4;
    dMwcr_r = if_else(K < K3, (Mwcr_true - Mwcr_r) / tau, 0);

    k7_updated = k7 * exp(-pA*max(1/Vf-1/Vfcr,0)) *(Mwcr_r/(Mg_A* (Nmu2+Nlambda2)/(Nmu1+Nlambda1)))^1.75;
    % k7_updated = k7 * exp(-pA*max(1/Vf-1/Vfcr,0)) *((6.25e+06/((Mg_A*(Nmu2+Nlambda2)/(Nmu1+Nlambda1))))^1.75);
    k7 = if_else(K >= K3 , k7_updated, k7);
end


% Differential equations
%dlambda0=k2*C_A*C_R-(k5*C_C+k6*C_E)*lambda0-k7*lambda0*lambda0;

dlambda0=k2*C_A*C_R-(k5*C_C+k6*C_E)*lambda0-k7*lambda0*lambda0;
dlambda1=k2*C_A*C_R+(k3+k4)*C_A*lambda0-(k4*C_A+k5*C_C+k6*C_E)*lambda1-k7*lambda0*lambda1;
dlambda2=k2*C_A*C_R+(k3+k4)*C_A*lambda0+2*k3*C_A*lambda1-(k4*C_A+k5*C_C+k6*C_E)*lambda2-k7*lambda0*lambda2;
dmu0=(k4*C_A+k5*C_C+k6*C_E)*lambda0+0.5*k7*lambda0*lambda0;
dmu1=(k4*C_A+k5*C_C+k6*C_E)*lambda1+0.5*k7*(lambda0*lambda1+lambda1*lambda0);
dmu2=(k4*C_A+k5*C_C+k6*C_E)*lambda2+0.5*k7*(lambda0*lambda2+2*lambda1*lambda1+lambda2*lambda0);
dC_R=2*f*k1*C_B-k2*C_A*C_R+(k5*C_C+k6*C_E)*lambda0;
dC_A=(1000/Mg_A)*Fgl/Vl-k2*C_A*C_R-(k3+k4)*lambda0*C_A;
dC_B=(1000/Mg_B)*FB/Vl-2*f*k1*C_B;
dC_C=-k5*lambda0*C_C;
dC_E=-k6*lambda0*C_E;

dNlambda0=dlambda0*Vl;
dNlambda1=dlambda1*Vl;
dNlambda2=dlambda2*Vl;
dNmu0=dmu0*Vl;
dNmu1=dmu1*Vl;
dNmu2=dmu2*Vl;
dN_R=dC_R*Vl;
dN_A=dC_A*Vl;
dN_B=dC_B*Vl;
dN_C=dC_C*Vl;
dN_E=dC_E*Vl;

dTr=(FA*Cp_A*(TA-Tr)+UA*(Tjacket-Tr)+deltaH*(k2*C_A*C_R+(k3+k4)*lambda0*C_A)*Vl-a*(Tr-Tamb)^b+Qstir)/(M_t*Cp + (Nmu1+Nlambda1)*(Mg_A/1000)*Cp_p+ NgA0*(Mg_A/1000)*Cp_A); % K/s
dNgA=(1000/Mg_A)*(FA-Fgl);
dP=R/(Vg/1000)*(Tr*dNgA+NgA*dTr); % kg/(m s^2)
dTjacket=(Fjacket*Cp*(Tcold-Tjacket)+Fjacket_h*Cp*(Thot-Tjacket)-UA*(Tjacket-Tr))/(Mj*Cp);
dMAin=FA;
dMBin=FB;
dMA_l=Fgl;
dVl=Fgl*0.06667;
% dTab=dTr+deltaH/(Cp*M_t)*dC_A-(dMAin+dMBin)*C_A*deltaH/(M_t^2*Cp);

% Only measurement
RA=(k2*C_A*C_R+(k3+k4)*lambda0*C_A)*3600*(Mg_A/1000)*Vl;%kg/h reaction rate of A
Mn=Mg_A*(mu1+lambda1)/(lambda0+mu0);%number-average molecular weight
Mw=Mg_A*(mu2+lambda2)/(mu1+lambda1);%weight-average molecular weight
PDI=Mw/max(Mn,1); % PDI

% lognormal distribution
sigma=sqrt(log(PDI));
mu=log(Mw./exp(sigma.^2./2));

% Cumulative Folry distribution  quasi-stationary
% tau = (k3*C_A)/(k3*C_A+k4*C_A+k5*C_C+k6*C_E+k7*lambda0);%(k3*C_A)/(k3*C_A+k4*C_A+k5*C_C+k6*C_E+k7);
% n=MX.sym('n',1);
% ff=-(1+(n+1)*(1-tau))*tau^(n+1)+(1+n*(1-tau))*tau^n;
% dWn=1/2*dmu1*ff*par.tf;

diff = vertcat(dNlambda0,dNlambda1,dNlambda2,dNmu0,dNmu1,dNmu2,dN_R,dN_A,dN_B,dN_C,dN_E,dTr,dNgA,dP,dTjacket,dMAin,dMBin,dMA_l,dVl);
x_var = vertcat(Nlambda0,Nlambda1,Nlambda2,Nmu0,Nmu1,Nmu2,N_R,N_A,N_B,N_C,N_E,Tr,NgA,P,Tjacket,MAin,MBin,MA_l,Vl);

if par.isFault==1
    diff = vertcat(dNlambda0,dNlambda1,dNlambda2,dNmu0,dNmu1,dNmu2,dN_R,dN_A,dN_B,dN_C,dN_E,dTr,dNgA,dP,dTjacket,dMAin,dMBin,dMA_l,dVl,dMwcr_r);
    x_var = vertcat(Nlambda0,Nlambda1,Nlambda2,Nmu0,Nmu1,Nmu2,N_R,N_A,N_B,N_C,N_E,Tr,NgA,P,Tjacket,MAin,MBin,MA_l,Vl,Mwcr_r);
    
end

d_var = vertcat(deltaH,UA,Thot,Tcold);
% u_var = vertcat(FA,FB,Fjacket,Fjacket_h);
u_var = vertcat(FA,FB,alpha);
k70_effective = k7 / exp(-E70/(R*Tr));  % Pre-exponential factor (shows Gel Effect modification)
y_var = vertcat(C_R,C_A,C_B,C_C,C_E,Tr,NgA,P,Tjacket,RA,Mn,Mw,PDI,MAin,MBin,Vf,k7*lambda0*lambda0,X,sigma,mu,k70_effective,Fgl);
% y_var = vertcat(Tr,P,RA,MAin,FA,FB,Fjacket);

% L = -10*(lambda0+mu0)*Vl+1e0*(Tr-Trsp).^2+1e2*FA^2+1e12*FB^2+1e-2*Fjacket^2+1e-2*Fjacket_h^2+(P-Psp).^2/1e4;

L = -20*dNmu1+2e3*(Tr-Trsp).^2+1e-1*(P-Psp).^2/1e4;%
par.istf=0;

Ltf=-20*Nmu1;
par.isg=0;


sys.f = Function('f',{x_var,u_var,d_var},{diff,L,Ltf},{'x','u','d'},{'xdot','qj','Ltf'});
sys.m = Function('m',{x_var,u_var,d_var},{y_var},{'x','u','d'},{'y'});


% sys.MWD = Function('MWD',{n,x_var,u_var,d_var},{dWn},{'n','x','u','d'},{'MWD'});

sys.ode = struct('x',x_var,'p',vertcat(u_var,d_var),'ode',diff,'quad',L);

% create CVODES integrator
sys.F = integrator('F','idas',sys.ode,struct('tf',par.tf));


sys.x = x_var;
sys.y = y_var;
sys.u = u_var;
sys.d = d_var;
sys.dx = diff;
sys.L = L;
state_index = struct( ...
    'Nlambda0', 1, ...
    'Nlambda1', 2, ...
    'Nlambda2', 3, ...
    'Nmu0', 4, ...
    'Nmu1', 5, ...
    'Nmu2', 6, ...
    'N_R', 7, ...
    'N_A', 8, ...
    'N_B', 9, ...
    'N_C', 10, ...
    'N_E', 11, ...
    'Tr', 12, ...
    'NgA', 13, ...
    'P', 14, ...
    'Tjacket', 15, ...
    'MAin', 16, ...
    'MBin', 17, ...
    'MA_l', 18, ...
    'Vl', 19);
if par.isFault==1
    state_index.Mwcr_r = 20;
end
input_index = struct('FA', 1, 'FB', 2, 'alpha', 3);
disturbance_index = struct('deltaH', 1, 'UA', 2, 'Thot', 3, 'Tcold', 4);
output_index = struct( ...
    'C_R', 1, ...
    'C_A', 2, ...
    'C_B', 3, ...
    'C_C', 4, ...
    'C_E', 5, ...
    'Tr', 6, ...
    'NgA', 7, ...
    'P', 8, ...
    'Tjacket', 9, ...
    'RA', 10, ...
    'Mn', 11, ...
    'Mw', 12, ...
    'PDI', 13, ...
    'MAin', 14, ...
    'MBin', 15, ...
    'Vf', 16, ...
    'kt', 17, ...
    'X', 18, ...
    'sigma', 19, ...
    'mu', 20, ...
    'k70_effective', 21, ...
    'Fgl', 22);
sys.audit = struct( ...
    'state_index', state_index, ...
    'input_index', input_index, ...
    'disturbance_index', disturbance_index, ...
    'output_index', output_index, ...
    'constants', struct( ...
        'Fmax_jc', Fmax_jc, ...
        'Cp_w', Cp, ...
        'M_liq', M, ...
        'M_reactor', M, ...
        'Cp_mix', Cp, ...
        'Cp_poly', Cp_p, ...
        'Cp_A', Cp_A, ...
        'Mg_A', Mg_A, ...
        'Mg_B', Mg_B, ...
        'NgA0', NgA0, ...
        'R_gas', R, ...
        'k10', k10, ...
        'k30', k30, ...
        'k40', k40, ...
        'k70', k70, ...
        'E10', E10, ...
        'E30', E30, ...
        'E40', E40, ...
        'E70', E70, ...
        'f_initiator', f, ...
        'Tcold_default', 283, ...
        'alpha_switch', alpha_switch));
% sys.nlcon = [nlcon];

%(lambda0,lambda1,lambda2,mu0, mu1, mu2, C_R, C_A, C_B, C_C, C_E, Tr,   NgA,  P,     Tjacket,  MAin,   MBin);
par.lbx = [0;      0;      0;       0;   0;   0;   0;   0;   0;   0;   0;   320;    0;  1.5e5;  283;      0;      0;    0;   2000;];
par.ubx = [inf;    inf;    inf;     inf; inf; inf; inf; inf; inf; inf; inf; 380;  10000; 1.6e6;  363;      Mmax_A; Mlim_B;  inf;  6000;];
%kg Feeding APS mass limit for changing the B feed rate
% FA,FB,Fj,Tjin
% par.lbu = [FA_min;FB_min;Fmin_jc;Fmin_jh;];
% par.ubu = [FA_max;FB_max;Fmax_jc;Fmax_jh;];
par.lbu = [FA_min;FB_min;-1];
par.ubu = [FA_max;FB_max; 1];

par.lbnlcon=0;
par.ubnlcon=0;


%(lambda0,lambda1,lambda2,mu0, mu1, mu2, C_R, C_A,    C_B,     C_C,   C_E,     Tr,   NgA,  P,     Tjacket, MAin,MBin);
par.x0 = [0;      0;      0;      0;   0;   0;   0;   C_A0*Vl0;   9.08E-08*Vl0;0.0021*Vl0;0.000252*Vl0;340;  NgA0;   P0; 340;    0;   0;     0; Vl0;];%340+172700/Cp*C_A0/M
% par.u0 = [FA_max;FB_min;Fmin_jc;Fmin_jh];
par.u0 = [FA_max;FB_min;0];
par.d0 = [272700;26870;363;283];
par.d0 = [272700*1.2;26870*0.8;363;283];

par.nx = numel(sys.x);
par.nu = numel(sys.u);
par.nd = numel(sys.d);


x_norm= [1.01391520272628e-08*4000
    0.00641597157351468*4000
    8191.96525456451*4000
    8.55985161634577e-06*4000
    5.69377231392881*4000
    7586277.75245485*4000
    1.30913280550926e-06*4000
    2.28534037943788*4000
    3.49851839611544e-07*4000
    0.00210000000000000*4000
    0.000252000000000000*4000
    351.858942899905
    8902.43968886215
    1602711.91550743
    347.635306357192
    3250.21793613436
    0.00168620000000030
    3000
    4000];

% u_norm=[0.694;1e-6;30;1;];
u_norm=[0.694;1e-6;1];

par.x_norm=x_norm;
par.u_norm=u_norm;


par.t0=[3975;10000;14000];
par.nt = numel(par.t0);

par.istf=0;


if par.isFault==1
    par.x0 = [0;      0;      0;      0;   0;   0;   0;   C_A0*Vl0;   9.08E-08*Vl0;0.0021*Vl0;0.000252*Vl0;340;  NgA0;   P0; 340;    0;   0;     0; Vl0; 1e6];
par.lbx = [0;      0;      0;       0;   0;   0;   0;   0;   0;   0;   0;   320;    0;  1.5e5;  283;      0;      0;    0;   2000; 1 ];
par.ubx = [inf;    inf;    inf;     inf; inf; inf; inf; inf; inf; inf; inf; 380;  10000; 1.6e6;  363;      Mmax_A; Mlim_B;  inf;  6000; inf];
x_norm= [1.01391520272628e-08*4000
    0.00641597157351468*4000
    8191.96525456451*4000
    8.55985161634577e-06*4000
    5.69377231392881*4000
    7586277.75245485*4000
    1.30913280550926e-06*4000
    2.28534037943788*4000
    3.49851839611544e-07*4000
    0.00210000000000000*4000
    0.000252000000000000*4000
    351.858942899905
    8902.43968886215
    1602711.91550743
    347.635306357192
    3250.21793613436
    0.00168620000000030
    3000
    4000
    6.26e6];
par.x_norm=x_norm;
end
