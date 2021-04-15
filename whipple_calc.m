%{ 
    How to use whipple_calc.m
    1) Edit your design projectile parameters
    2) Edit your Wall/Bumper spacing and material (density/rear wall yield stress)
    3) Edit your Nextel and Kevlar material (data provided in script)
    4) Profit
    This script can spit out:
    - Number of Nextel/Kevlar layers, rounded up for redundency
    - Wall and bumper thickness in cm
    - Critical projectile diameter in cm (vs velocity)

    This line is to demonstrate that git works
%}
clc
% Projectile parameters
theta=0; % impact angle (degree)
% p_flux=0.05;% #per m2 per year
% d_design=interp1(a400(:,2),a400(:,1),p_flux)*100; % Design projectile diameter (cm)
d_design=.05; 
rho_p= 2.8; % projectille density (g/cm3)
v=22; % projectile velocity (km/s)
v_n=v*cosd(theta); % projectile velocity normal to bumper(km/s)

% Wall/Bumper parameters
s=10; %overall spacing (cm)
rho_b= 2.71; % bumper density (g/cm3) Al-6061-T6
rho_w=2.84; % rear wall density (g/cm3) Al-2219-T87
sig= 0.145038*390; % rear wall yield stress (ksi) 2219-T87

% Stuffing parameters
%%% Nextel: 0.1 (AF62), 0.027 (AF10)
%%% Kevlar: 0.023 (KM2-705), 0.032 (22-FDI129 aka 710)
c_nk=0.23; % CONSTANT, DO NOT CHANGE
m_kev=0.032; m_nex=0.1;

%%%% Regular Whipple Shielding
% d_crit=whip_crit(t_w,rho_b,rho_p,v_n,s,sig);
% [t_b,t_w]=whip_sizing(rho_b,rho_p,v_n,s,sig,d);

% Areal density and number of Nextel/Kevlar layers
[m_nk,no_nex,no_kev]=nexkev_sizing(c_nk,d_design,rho_p,m_kev,m_nex);
% m_nk=m_nk+m_kev;
% Critical wall/bumper sizing for design projectile
[t_b,t_w]=stuff_sizing(rho_w,rho_b,rho_p,v_n,s,sig,d_design,c_nk,theta)

%%%% Test data for t_w=0.48cm, t_wb.2cm, 6 Nex AF62 + 6 Kev 710, S=10.7cm
%%%% d_crit=stuff_crit(.48,.2,2.84,2.713,2.796,6,10.7,sig,theta,0.796)

% Critical projectile diameter for design wall/bumper
% d_crit=stuff_crit(t_w,t_b,rho_w,rho_b,rho_p,2,s,sig,theta,m_nk)
[d_crit,mshield]=stuff_crit(round(t_w,2),round(t_b,2),rho_w,rho_b,rho_p,22,s,sig,theta,m_nk)

vv=1:.1:15;
for i=1:length(vv)
    [dd(i),m_shield(i)]=stuff_crit(0.24,0.05,rho_w,rho_b,rho_p,vv(i),s,sig,theta,m_nk);
%     dd(i)=stuff_crit(.48,.2,2.84,2.713,2.796,vv(i),10.7,sig,theta,0.796);
%     dd(i)=stuff_crit(0.25,0.1,rho_w,rho_b,rho_p,2,s,sig,theta,m_nk)
end
plot(vv,dd)
ylim([0,2])
%% Whippple functions
function d_crit=whip_crit(t_w,rho_b,rho_p,v_n,s,sig)
if v_n>=3+7.7
    d_crit=3.918*t_w^(2/3)*rho_p^(-1/3)*rho_b^(-1/9)*v_n^(-2/3)*s^(1/3)*(sig/70)^(1/3);
else
    d_crit=((t_w*(sig/40)^0.5+t_b)/(0.6*rho_p^0.5*v_n^2/3))^(18/19);
end
end
function [t_b,t_w]=whip_sizing(rho_b,rho_p,v_n,s,sig,d)
c_w=0.16;
M_p=4/3.*(d/2).^3*pi*rho_p;
if s/d>=30
    c_b=0.2;
else
    c_b=0.25;
end
t_b=c_b.*d*rho_p/rho_b;
t_w=c_w.*d.^0.5*(rho_p*rho_b)^(1/6).*(M_p).^(1/3)*v_n/s^0.5*(70/sig)^0.5;
end
%% Stuffed Whipple functions
function [m_nk,no_nex,no_kev]=nexkev_sizing(c_nk,d,rho_p,m_kev,m_nex)
m_nk_design=c_nk*d*rho_p;
no_nex=layer_round(m_nk_design*0.75/m_nex);
no_kev=layer_round(m_nk_design*0.25/m_kev);
m_nk=no_nex*m_nex+no_kev*m_kev;
end
function no_layer=layer_round(inpt)
if (inpt-round(inpt))>=0.15
    no_layer=ceil(inpt);
else
    no_layer=round(inpt);
end
end
function [t_b,t_w]=stuff_sizing(rho_w,rho_b,rho_p,v_n,s,sig,d,theta,m_nk)
c_w=8.84;
c_o=0.38; c_b=0.15;
M_p=4/3.*(d/2).^3*pi*rho_p;
t_b=c_b.*d*rho_p/rho_b;
t_w=c_w*(c_o*d*rho_p/(t_b*rho_p+m_nk))^1.1*M_p^(1/3)*v_n*cosd(theta)^0.5*1/rho_w*s^-2*(sig/40)^(-.5);
end
function [d_crit,m_shield]=stuff_crit(t_w,t_b,rho_w,rho_b,rho_p,v_n,s,sig,theta,m_nk)
m_b=rho_b*t_b; m_w=rho_w*t_w;
m_shield=m_b+m_nk+m_w;
m_b_total=m_b+m_nk;
v=v_n/cosd(theta);
cL=0.37;
if v>=6.5*cosd(theta)^(-0.75) % High vel impact
    if m_nk>=0.25*m_shield && m_nk<=0.35*m_shield
        k_h_sw=0.6;
    elseif m_nk>=0.1*m_shield && m_nk<=0.15*m_shield
        k_h_sw=0.45;
    else
        k_h_sw=0.53;
    end
%     fprintf('High velocity impact\n')
    d_crit=k_h_sw*(t_w*rho_w)^(1/3)*rho_p^(-1/3)*(sig/40)^(1/6)*v^(-1/3)*cosd(theta)^(-0.5)*s^(2/3);
elseif v<=2.6*cosd(theta)^(-0.5)
    k_l_sw=2.35;
%     fprintf('Low Velocity impact\n')
    d_crit=k_l_sw*v^(-2/3)*cosd(theta)^(-4/3)*rho_p^(-0.5)*(t_w*(sig/40)^0.5+cL*m_b_total);
else
    if m_nk>=0.25*m_shield && m_nk<=0.35*m_shield
        k_h_sw=0.321;
    elseif m_nk>=0.1*m_shield && m_nk<=0.15*m_shield
        k_h_sw=0.241;
    else
        k_h_sw=0.281;
    end
    k_l_sw=1.243;
%     fprintf('Intermediate velocity impact\n')
    d_crit=(k_l_sw*(t_w*(sig/40)^0.5+cL*m_b_total)/(cosd(theta)*rho_p^0.5))*...
        (6.5*cosd(theta)^(-0.75)-v)/(6.5*cosd(theta)^(-0.75)-2.6*cosd(theta)^(-0.5))...
        +...
        k_h_sw*(t_w*rho_w)^(1/3)*rho_p^(-1/3)*(sig/40)^(1/6)*cosd(theta)^(-0.25)*s^(2/3)*...
        (v-2.6*cosd(theta)^(-0.5))/(6.5*cosd(theta)^(-0.75)-2.6*cosd(theta)^(-0.5));        
end
end