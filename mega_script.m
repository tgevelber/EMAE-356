%%% insert your code here
%% MMOD Shielding model, Single layer
% Aluminum
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