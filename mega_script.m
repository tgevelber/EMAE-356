%%% insert your code here
%%% MMOD Shielding model, Single layer
%% Aluminum
d= % Projectile diameter (cm)
rho_p= % projectille density (g/cm3)
rho_t= % target density (g/cm3)
bhn= % target Brinell hardness
c_t= % Speed of sound of target (km/s)
theta= % impact angle from target normal (deg)
v= % projectile velocity

% Density ratio coefficient
if rho_p/rho_t>=1.5
    n=2/3;
else
    n=0.5;
end
% penetration depth (cm)
p_inf=5.24*d^(19/18)*bnh^(-0.25)*(rho_p/rho_t)^n*(v*cosd(theta)/c_t)^(2/3);

% thickness limit (cm) for incipient, detached, perforation
k=[3 2.2 1.8];
t=p_inf*k;

% projectile max diameter for specific wall config
d_c=((t*bhn^0.25*(rho_t/rho_p)^0.5)/...
    (k*5.24*(v*cosd(theta)/c_t)^2/3))^(18/19);

%% Whipple (double Al)
s= 11;%overall spacing (cm)
d= 1;%projectile diameter (cm)
if s/d>=30
    c_b=0.2;
else
    c_b=0.25;
end
c_w=0.16;
sig= 0.145038*71;% rear wall yield stress
rho_p= 2.8;% projectille density (g/cm3)
rho_b= 2.74;% bumper density (g/cm3)
v_n=7;
% t_b=c_b*m_p/rho_b;
M_p=4/3*(d/2)^3*pi*rho_p;

% Bumper thickness(cm)
t_b=c_b*d*rho_p/rho_b; 
% Rear wall thickness(cm)
t_w=c_w*d^0.5*(rho_p*rho_b)^(1/6)*(M_p)^(1/3)*v_n/s^0.5*(70/sig)^0.5;