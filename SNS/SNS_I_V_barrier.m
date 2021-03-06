% SNS Ballistic Device
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
t = 1;

% transmission paramter
% t value for the coupling element
trans = 1.0

% mu : Device Fermi Level
mu = 2;

% kT : k_B * Temperature(K)
% kT units : eV
% kT default value : 0.026 eV (Room temperature)
kT = 0.00001;

% phase difference betweem the two superconductors
phi = 0;

% Superconducting order paramter
Delta1 = 0;
Delta2 = 1e-2*exp(1j*phi);

% N_V : number of points in the voltage grid
N_V = 100;
V_vec = 3e-2 * linspace(0,1,N_V);
I_vec = zeros(4,length(V_vec));

for ii = 1:length(V_vec)
    V = V_vec(ii)
    
    % mu1 : electrochemical potential of the left contact
    mu1 = mu + V;
    % mu2 : electrochemical potential of the right contact
    mu2 = mu;

    U = 0;
    I_E = @(E) calculate_I_E(E,t,trans,U,mu,mu1,mu2,Delta1,Delta2,kT,eta,N_D);
    I_vec(1,ii) = integral(I_E,-1.5*abs(V),1.5*abs(V),'AbsTol',1e-6,'ArrayValued',true);
    %U = 0.5;
    %I_E = @(E) calculate_I_E(E,t,trans,U,mu,mu1,mu2,Delta1,Delta2,kT,eta,N_D);
    %I_vec(2,ii) = integral(I_E,-1.5*abs(V),1.5*abs(V),'AbsTol',1e-6,'ArrayValued',true);
    %U = 1;
    %I_E = @(E) calculate_I_E(E,t,trans,U,mu,mu1,mu2,Delta1,Delta2,kT,eta,N_D);
    %I_vec(3,ii) = integral(I_E,-1.5*abs(V),1.5*abs(V),'AbsTol',1e-6,'ArrayValued',true);
    %U = 2;
    %I_E = @(E) calculate_I_E(E,t,trans,U,mu,mu1,mu2,Delta1,Delta2,kT,eta,N_D);
    %I_vec(4,ii) = integral(I_E,-1.5*abs(V),1.5*abs(V),'AbsTol',1e-6,'ArrayValued',true);
    
end
dV = V_vec(2) - V_vec(1);
figure(1)
plot(V_vec(2:length(V_vec))./Delta2,diff(I_vec(1,:))./dV,'linewidth',2.5);
%hold on;
%plot(V_vec(2:length(V_vec))./Delta2,diff(I_vec(2,:))./dV,'linewidth',2.5);
%plot(V_vec(2:length(V_vec))./Delta2,diff(I_vec(3,:))./dV,'linewidth',2.5);
%plot(V_vec(2:length(V_vec))./Delta2,diff(I_vec(4,:))./dV,'linewidth',2.5);
%hold off;
set(gca,'FontSize',20)
ylim([0,5])
xlabel('$\frac{eV}{\Delta}$','interpreter','latex','fontsize',24);
ylabel('Conductance G ($\frac{e^2}{h}$)','interpreter','latex','fontsize',24);
title(['G(V)'],'interpreter','latex','fontsize',20);
%legend('U = 0','U = 0.5','U = 1.0','U = 2.0');
