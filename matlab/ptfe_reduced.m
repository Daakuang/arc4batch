function [sys,par] = ptfe_reduced(par)

% industrial batch polymerization reactor PTFE
%
% Written by: Dr. Chenchen Zhou, July. 2024 ZJU

%u=[FA,FB,Fj,Tjin] manipulated variable
%y=[C_R,C_A,C_B,C_C,C_E,Tr,NgA,P,Tj,RA,Mn,Mw,PDI] Measured variable

import casadi.*

Nlambda0=MX.sym('Nlambda0',1); %  moment of active polymer
% Nlambda1=MX.sym('Nlambda1',1); %  moment of active polymer
% Nlambda2=MX.sym('Nlambda2',1); %  moment of active polymer
% Nmu0=MX.sym('Nmu0',1); %  moment of dead polymer
% Nmu1=MX.sym('Nmu1',1); %  moment of dead polymer
% Nmu2=MX.sym('Nmu2',1); %  moment of dead polymer
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

Fjacket=(Fmax_jc-Fmin_jc)*max(-alpha,0); % kg/s Flow rate of jacket water
Fjacket_h=(Fmax_jh-Fmin_jh)*max(alpha,0); % kg/s Flow rate of jacket water

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

k10=1.13e16*1e1; %1/s Pre-reference factor
k20=3.62e13*1e2;
k30=5.49e6*1e1;1.249e7;
k40=3.32e6*3e1;
k50=1051*5e2*6;
k60=1051*1e3*0;
k70=3.38e4*1e5;1.38e4;
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
% Mmax_A=3250; %kg Maximum feed mass of A

Trsp=351; %K Reactor temperature set point
Psp=1500000; %Pa Reactor gas pressure set point
% Pact_Tc=1000000; %Pa Pressure standard for starting the temperature control

% algebraic equations

lambda0=Nlambda0/Vl; %  moment of active polymer
% lambda1=Nlambda1/Vl; %  moment of active polymer
% lambda2=Nlambda2/Vl; %  moment of active polymer
% mu0=Nmu0/Vl; %  moment of dead polymer
% mu1=Nmu1/Vl; %  moment of dead polymer
% mu2=Nmu2/Vl; %  moment of dead polymer
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
M_t  = M + NgA0*(Mg_A/1000) + MAin + MBin;
% Vl=Vl0+MA_l*0.6667;
Vg=V-Vl;
f = 1;
% Glass effect
% k3=k3*exp(-pB*0.5*((1/Vf-1/Vfcr2)+sqrt((1/Vf-1/Vfcr2)^2+1e-5)));
% k3=k3*exp(-pB*max(1/Vf-1/Vfcr2,0));

% Cage effect
% f = f*exp(-pC*0.5*((1/Vf-1/Vfcr3)+sqrt((1/Vf-1/Vfcr3)^2+1e-5)));
% f = f0*exp(-pC*max(1/Vf-1/Vfcr3,0));
% Gel effect
% K = Mw^m*exp(pA/Vf)-K3;>0
if par.isFault==1
    % M_t  = M + NgA0*(Mg_A/1000) + MAin + MBin;
    % Vl=Vl0+MA_l*0.6667;
    % Fix for missing MA_l state: Derive from Vl (Vl = Vl0 + MA_l * 0.06667)
    MA_l = (Vl - Vl0) / 0.06667; % Approximate relation based on dVl=0.06667*Fgl
    X=(C_A0*Vl0+MA_l*(1000/Mg_A)-C_A*Vl)/(C_A0*Vl+MA_l*(1000/Mg_A));
    % For 3 effects
    % Mwcr=6e7;
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
    ap = 4.8e-4;
    aM = 1e-3;

    % free-volume
    phi_p=(1+epsilon)*X/(1+epsilon*X);
    phi_M=1-phi_p;
    Vf = 0.025+ap*(Tr-Tgp)*phi_p+aM*(Tr-TgM)*phi_M;
    % k5=k5*exp(-pA*0.5*((1/Vf-1/Vfcr)+sqrt((1/Vf-1/Vfcr)^2+1e-5)));
    % k6=k6*exp(-pA*0.5*((1/Vf-1/Vfcr)+sqrt((1/Vf-1/Vfcr)^2+1e-5)));
    k7=k7*exp(-pA*0.5*((1/Vf-1/Vfcr)+sqrt((1/Vf-1/Vfcr)^2+1e-5)));
end


% Differential equations
%dlambda0=k2*C_A*C_R-(k5*C_C+k6*C_E)*lambda0-k7*lambda0*lambda0;

dlambda0=k2*C_A*C_R-(k5*C_C)*lambda0-k7*lambda0*lambda0;
% dlambda1=k2*C_A*C_R+(k3+k4)*C_A*lambda0-(k4*C_A+k5*C_C+k6*C_E)*lambda1-k7*lambda0*lambda1;
% dlambda2=k2*C_A*C_R+(k3+k4)*C_A*lambda0+2*k3*C_A*lambda1-(k4*C_A+k5*C_C+k6*C_E)*lambda2-k7*lambda0*lambda2;
% dmu0=(k4*C_A+k5*C_C+k6*C_E)*lambda0+0.5*k7*lambda0*lambda0;
% dmu1=(k4*C_A+k5*C_C+k6*C_E)*lambda1+0.5*k7*(lambda0*lambda1+lambda1*lambda0);
% dmu2=(k4*C_A+k5*C_C+k6*C_E)*lambda2+0.5*k7*(lambda0*lambda2+2*lambda1*lambda1+lambda2*lambda0);
dC_R=2*f*k1*C_B-k2*C_A*C_R+(k5*C_C)*lambda0;
dC_A=(1000/Mg_A)*Fgl/Vl-k2*C_A*C_R-(k3+k4)*lambda0*C_A;
dC_B=(1000/Mg_B)*FB/Vl-2*f*k1*C_B;
dC_C=-k5*lambda0*C_C;

dNlambda0=dlambda0*Vl;
% dNlambda1=dlambda1*Vl;
% dNlambda2=dlambda2*Vl;
% dNmu0=dmu0*Vl;
% dNmu1=dmu1*Vl;
% dNmu2=dmu2*Vl;
dN_R=dC_R*Vl;
dN_A=dC_A*Vl;
dN_B=dC_B*Vl;
dN_C=dC_C*Vl;

dTr=(FA*Cp_A*(TA-Tr)+UA*(Tjacket-Tr)+deltaH*(k2*C_A*C_R+(k3+k4)*lambda0*C_A)*Vl-a*(Tr-Tamb)^b+Qstir)/(M*Cp+(NgA0*(Mg_A/1000) + MAin)*Cp_A);
dNgA=(1000/Mg_A)*(FA-Fgl);
dP=R/(Vg/1000)*(Tr*dNgA+NgA*dTr); % kg/(m s^2)
dTjacket=(Fjacket*Cp*(Tcold-Tjacket)+Fjacket_h*Cp*(Thot-Tjacket)-UA*(Tjacket-Tr))/(Mj*Cp);
dMAin=FA;
dMBin=FB;
% dMA_l=Fgl;
dVl=Fgl*0.06667;
% dTab=dTr+deltaH/(Cp*M_t)*dC_A-(dMAin+dMBin)*C_A*deltaH/(M_t^2*Cp);

% Only measurement
RA=(k2*C_A*C_R+(k3+k4)*lambda0*C_A)*Vl;%kg/h reaction rate of A *3600*(Mg_A/1000)
% Mn=Mg_A*(mu1+lambda1)/(lambda0+mu0);%number-average molecular weight
% Mw=Mg_A*(mu2+lambda2)/(mu1+lambda1);%weight-average molecular weight
% PDI=Mw/max(Mn,1); % PDI

% lognormal distribution
% sigma=sqrt(log(PDI));
% mu=log(Mw./exp(sigma.^2./2));

% Cumulative Folry distribution  quasi-stationary
% tau = (k3*C_A)/(k3*C_A+k4*C_A+k5*C_C+k6*C_E+k7*lambda0);%(k3*C_A)/(k3*C_A+k4*C_A+k5*C_C+k6*C_E+k7);
% n=MX.sym('n',1);
% ff=-(1+(n+1)*(1-tau))*tau^(n+1)+(1+n*(1-tau))*tau^n;
% dWn=1/2*dmu1*ff*par.tf;



diff = vertcat(dNlambda0,dN_R,dN_A,dN_B,dN_C,dTr,dNgA,dP,dTjacket,dMAin,dMBin,dVl);
%              (Nlambda0,Nlambda1,Nlambda2,Nmu0,Nmu1,Nmu2,N_R,N_A,N_B,N_C,N_E,Tr,NgA,P,Tjacket,MAin,MBin,MA_l,Vl)
x_var = vertcat(Nlambda0,                                 N_R,N_A,N_B,N_C,    Tr,NgA,P,Tjacket,MAin,MBin,     Vl);
d_var = vertcat(deltaH,UA,Thot,Tcold);
% u_var = vertcat(FA,FB,Fjacket,Fjacket_h);
u_var = vertcat(FA,FB,alpha);
y_var = vertcat(Tr,P,Tjacket);
% y_var = vertcat(Tr,P,RA,MAin,FA,FB,Fjacket);

% L = -10*(lambda0+mu0)*Vl+1e0*(Tr-Trsp).^2+1e2*FA^2+1e12*FB^2+1e-2*Fjacket^2+1e-2*Fjacket_h^2+(P-Psp).^2/1e4;

% L = -20*RA+2e3*(Tr-Trsp).^2+1e-1*(P-Psp).^2/1e4;%
% t=0.5*(Mmax_A-MAin+sqrt((Mmax_A-MAin).^2+1e-5));
% t=tanh((Mmax_A-MAin));
if par.SMS==1
    t=1/(1+exp(-10*(Mmax_A-MAin)+10));
else
    t=MAin<Mmax_A;
end
% L = (-1*RA+1e5/360*(P/Psp-1).^2)*t+1e2/360*(Tr-Trsp).^2;%
L = -1*RA+(3e5/360*(P/Psp-1).^2)*t+3e2/360*(Tr-Trsp).^2;%
% L = (-20*RA+1e-1*(P-Psp).^2/1e4)*(MAin<Mmax_A)+2e3*(Tr-Trsp).^2;%
%L = (-1*RA+1/40*1e-7*(P-Psp).^2)*t+2e2*(Tr-Trsp).^2;% N=300 t=1/(1+exp(-10*(Mmax_A-MAin)+10));
par.istf=0;

% Ltf=-20*Nmu1;
% t1=1./(1+exp(-2*(Trsp-Tr)+10));
par.isg=0;
% % nlcon = 2e2/360*(Tr-Trsp)*dTr*t1;
% par.lbnlcon=-inf;
% par.ubnlcon=0;

if ~par.isg
    sys.f = Function('f',{x_var,u_var,d_var},{diff,L},{'x','u','d'},{'xdot','qj'});
    sys.m = Function('m',{x_var,u_var,d_var},{y_var},{'x','u','d'},{'y'});
else
    sys.f = Function('f',{x_var,u_var,d_var},{diff,L,nlcon},{'x','u','d'},{'xdot','qj','gj'});
    sys.m = Function('m',{x_var,u_var,d_var},{y_var},{'x','u','d'},{'y'});
end


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
% sys.nlcon = [nlcon];

%(lambda0,lambda1,lambda2,mu0, mu1, mu2, C_R, C_A, C_B, C_C, C_E, Tr,   NgA,  P,     Tjacket,  MAin,   MBin);
par.lbx = [0;   0;   0;   0;   0;   320.15;    0;  1.5e5;  283;      0;      0;    2000;];
par.ubx = [inf; inf; inf; inf; inf; 373.15;  10000; 1.6e6;  373;      Mmax_A; Mlim_B;  6000;];
% par.ubx = [1.01391520272628e-08*8000;1.30913280550926e-06*8000;2.28534037943788*8000;3.49851839611544e-07*8000;0.00210000000000000*8000; 373.15;  10000; 1.6e6;  373;      Mmax_A; Mlim_B;  6000;];
%kg Feeding APS mass limit for changing the B feed rate
% FA,FB,Fj,Tjin
% par.lbu = [FA_min;FB_min;Fmin_jc;Fmin_jh;];
% par.ubu = [FA_max;FB_max;Fmax_jc;Fmax_jh;];
par.lbu = [FA_min;FB_min;-1];
par.ubu = [FA_max;FB_max; 1];

% par.lbnlcon=0;
% par.ubnlcon=0;


%(lambda0,lambda1,lambda2,mu0, mu1, mu2, C_R, C_A,    C_B,     C_C,   C_E,     Tr,   NgA,  P,     Tjacket, MAin,MBin);
par.x0 = [0;   0;   C_A0*Vl0;   9.08E-08*Vl0;0.0021*Vl0;340;  NgA0;   P0; 340;    0;   0;   Vl0;];%340+172700/Cp*C_A0/M
% par.u0 = [FA_max;FB_min;Fmin_jc;Fmin_jh];
par.u0 = [FA_max;FB_max;1];
par.d0 = [272700;26870;363;283];
par.d0 = [272700*1.2;26870*0.8;363;283];

par.nx = numel(sys.x);
par.nu = numel(sys.u);
par.nd = numel(sys.d);

%(Nlambda0,N_R,N_A,N_B,N_C,Tr,NgA,P,Tjacket,MAin,MBin,MA_l,Vl)
%(Nlambda0,Nlambda1,Nlambda2,Nmu0,Nmu1,Nmu2,N_R,N_A,N_B,N_C,N_E,Tr,NgA,P,Tjacket,MAin,MBin,MA_l,Vl)
% x_norm= [1.01391520272628e-08*4000
%     1.30913280550926e-06*4000
%     2.28534037943788*4000
%     3.49851839611544e-07*4000
%     0.00210000000000000*4000
%     351.858942899905
%     8902.43968886215
%     1602711.91550743
%     347.635306357192
%     3250.21793613436
%     0.00168620000000030
%     4000];
%

x_norm = [0.00015
    0.05
    9000
    0.0016
    10
    365
    1075
    1554000
    350
    3250
    0.002
    4225];

x_norm = [0.00015
    0.05
    9000
    0.0004
    10
    365
    1075
    1554000
    350
    3250
    0.002*1e-2
    4225];

x_norm1 = [0.00015
    0.05
    9000
    1e-7
    10
    365
    1075
    1554000
    350
    3250
    0.002
    4225];


% u_norm=[0.694;1e-6;30;1;];
u_norm=[FA_max;FB_max;1];

par.x_norm=x_norm;%./x_norm;
par.u_norm=u_norm;%./u_norm;
par.x_norm1=x_norm1;%./x_norm;


par.t0=[3975;10000;14000];
par.nt = numel(par.t0);

par.istf=0;
