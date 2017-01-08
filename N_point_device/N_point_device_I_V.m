% I-V Charactersitics for a N-point channel device
% Author: Sandesh Kalantre

% tight-binding parameter [eV]
t = 1.0;

% 0+ for the iteration to converge to one of the roots
eta = 1e-8;

% N_D : number of points in the device region
N_D = 100;

%Fermi level
mu = 2.0;

% Temperature kT : [eV]
kT = 0.001;

% N_V: number of points in the voltage grid
N_V = 100;
V_vec = 1*linspace(-1,1,N_V);
I_vec = zeros(1,length(V_vec));
for ii = 1:length(V_vec)
   V = V_vec(ii)
   
   mu1 = mu+V/2;
   mu2 = mu-V/2;
   
   I_E = @(E) calculate_I_E(E,t,mu1,mu2,kT,eta,N_D);
   I_vec(ii) = quadv(I_E,mu-2*abs(V),mu+2*abs(V),1e-5);
end

figure(1)
plot(V_vec,I_vec,'linewidth',1.5);
xlabel('Voltage (V)','fontsize',16,'interpreter','latex');
ylabel('Current (I)','fontsize',16,'interpreter','latex');
title('I-V for 100 point ballistic device','fontsize',16,'interpreter','latex');