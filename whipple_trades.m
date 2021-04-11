clean
%%% Wall thickness vs projectile diameter
s= 20;%overall spacing (cm)
d= 0.01:0.01:1;%projectile diameter (cm)
for i=1:length(d)
if s/d(i)>=30
    c_b(i)=0.2;
else
    c_b(i)=0.25;
end
end
c_w=0.16;
sig= 0.145038*71;% rear wall yield stress (ksi)2
rho_p= 2.8;% projectille density (g/cm3)
rho_b= 2.78;% bumper density (g/cm3)
v_n=22; % projectile speed normal to bumper(km/s)
% t_b=c_b*m_p/rho_b;
M_p=4/3.*(d/2).^3*pi*rho_p;

% Bumper thickness(cm)
t_b=c_b.*d*rho_p/rho_b; 
% Rear wall thickness(cm)
t_w=c_w.*d.^0.5*(rho_p*rho_b)^(1/6).*(M_p).^(1/3)*v_n/s^0.5*(70/sig)^0.5;

% critical particle sizing
% d_c=3.918*t_w^(2/3)*rho_p^(-1/3)*rho_b^(-1/9)*
figure (1)
loglog(d,t_b)
hold on
loglog(d,t_w)
title('Whipple sizing vs projectile diameter @22km/s, S=20cm')
xlabel('Projectile Diameter (cm)')
ylabel('Thickness (cm)')
legend('Bumper', 'Rear wall')
%% Wall thickness vs spacing for d=.5
clear
s= 10:.1:30;%overall spacing (cm)
d= .3;%projectile diameter (cm)
for i=1:length(s)
if s(i)/d>=30
    c_b(i)=0.2;
else
    c_b(i)=0.25;
end
end
c_w=0.16;
sig= 0.145038*71;% rear wall yield stress (ksi)2
rho_p= 2.8;% projectille density (g/cm3)
rho_b= 2.78;% bumper density (g/cm3)
v_n=22; % projectile speed normal to bumper(km/s)
% t_b=c_b*m_p/rho_b;
M_p=4/3.*(d/2).^3*pi*rho_p;

% Bumper thickness(cm)
t_b=c_b.*d*rho_p/rho_b; 
% Rear wall thickness(cm)
t_w=c_w.*d.^0.5*(rho_p*rho_b)^(1/6).*(M_p).^(1/3).*v_n./s.^0.5*(70/sig)^0.5;

% critical particle sizing
% d_c=3.918*t_w^(2/3)*rho_p^(-1/3)*rho_b^(-1/9)*
figure (2)
plot(s,t_b)
hold on
plot(s,t_w)
plot([20,20],[0,.8])
title('Whipple sizing vs spacing @0.3cm, 22km/s')
xlabel('Layer spacing (cm)')
ylabel('Thickness (cm)')
legend('Bumper', 'Rear wall')
%%
% %%% MMOD Shielding model, Single layer
% %%Aluminum
% % d= % Projectile diameter (cm)
% % rho_p= % projectille density (g/cm3)
% % rho_t= % target density (g/cm3)
% % bhn= % target Brinell hardness
% % c_t= % Speed of sound of target (km/s)
% % theta= % impact angle from target normal (deg)
% % v= % projectile velocity
% 
% % Density ratio coefficient
% if rho_p/rho_t>=1.5
%     n=2/3;
% else
%     n=0.5;
% end
% % penetration depth (cm)
% p_inf=5.24*d^(19/18)*bnh^(-0.25)*(rho_p/rho_t)^n*(v*cosd(theta)/c_t)^(2/3);
% 
% % thickness limit (cm) for incipient, detached, perforation
% k=[3 2.2 1.8];
% t=p_inf*k;
% 
% % projectile max diameter for specific wall config
% d_c=((t*bhn^0.25*(rho_t/rho_p)^0.5)/...
%     (k*5.24*(v*cosd(theta)/c_t)^2/3))^(18/19);

