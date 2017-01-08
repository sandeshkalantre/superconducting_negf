% NNN Ballistic Device
% Author: Sandesh Kalantre

% Calculation of current and conductance vs bias

% eta : 0+ for calculating retarded Green's functions
eta = 1e-8;

% Fundamental Constants
q = 1;
hbar = 1;
    
% number of points in the device (channel)
N_D = 2;

% t0 : tight binding parameter = hbar^2/(2 m a^2)
% t0 units : eV
t0 = 1;

% mu : Device Fermi Level
mu = 2.0;

% kT : k_B * Temperature(K)
% kT units : eV
% kT default value : 0.026 eV (Room temperature)
kT = 0.0001;

% N_V : number of points in the voltage grid
N_V = 100;
V_vec = 1 * linspace(-1,1,N_V);
I_vec = zeros(1,length(V_vec));

for ii = 1:length(V_vec)
    V = V_vec(ii)
    
    % mu1 : electrochemical potential of the left contact
    mu1 = mu + V/2;
    % mu2 : electrochemical potential of the right contact
    mu2 = mu - V/2;

    I_E = @(E) calculate_I_E(E,t,mu1,mu2,kT,eta,N_D);
    I_vec(ii) = quadv(I_E,-mu-abs(V/2),-mu+abs(V/2),1e-7) + quadv(I_E,mu-abs(V/2),mu+abs(V/2),1e-7);
end

figure(1)
plot(V_vec,I_vec,'linewidth',2.0);
set(gca,'FontSize',20)
xlabel('V','interpreter','latex','fontsize',20);
ylabel('Current (A)','interpreter','latex','fontsize',20);
title(['I-V for ' num2str(N_D) ' device point(s)'],'interpreter','latex','fontsize',20);
