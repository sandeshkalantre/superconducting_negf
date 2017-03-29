% SNS Ballistic Device
% Author: Sandesh Kalantre

% Calculation of current and conductance vs bias

% eta : 0+ for calculating retarded Green's functions
eta = 1e-6;

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
trans = 1;

% mu : Device Fermi Level
mu = 2;

% kT : k_B * Temperature(K)
% kT units : eV
% kT default value : 0.026 eV (Room temperature)
kT = 0.00001;

% phase difference betweem the two superconductors
phi = 0;

% Superconducting order paramter
Delta1 = 1e-2;
Delta2 = 1e-2*exp(1j*phi);

% N_V : number of points in the voltage grid
N_V = 100;
V_vec = 2e-2 * linspace(-1,1,N_V);
I_vec = zeros(4,length(V_vec));

for ii = 1:length(V_vec)
    V = V_vec(ii)
    
    % mu1 : electrochemical potential of the left contact
    mu1 = mu + V;
    % mu2 : electrochemical potential of the right contact
    mu2 = mu;

    U = 0;
    DOS_E = @(E) calculate_DOS(E,t,trans,U,mu,mu1,mu2,Delta1,Delta2,kT,eta,N_D);
    E_vec = 1.5*Delta1 * linspace(-1,1,100);
    DOS = arrayfun(DOS_E,E_vec);
    figure(2)
    plot(E_vec,DOS,'linewidth',2.5);
    drawnow;
    
    I_E = @(E) calculate_I_E(E,t,trans,U,mu,mu1,mu2,Delta1,Delta2,kT,eta,N_D);
    I_vec(1,ii) = integral(I_E,-1.5*abs(V),1.5*abs(V),'AbsTol',1e-8,'ArrayValued',true);
    
end
dV = V_vec(2) - V_vec(1);
%figure(1)
%conductance
%plot(V_vec(2:length(V_vec))./Delta2,diff(I_vec(1,:))./dV,'linewidth',2.5);
%set(gca,'FontSize',20)
%xlabel('$\frac{eV}{\Delta}$','interpreter','latex','fontsize',24);
%ylabel('Conductance G ($\frac{e^2}{h}$)','interpreter','latex','fontsize',24);
%title(['G(V)'],'interpreter','latex','fontsize',24);
%legend('U = 0','U = 0.5','U = 1.0','U = 2.0');

%IV
figure(1)
plot(V_vec./Delta1,I_vec,'linewidth',2.5);
set(gca,'FontSize',20)
xlabel('$\frac{eV}{\Delta_2}$','interpreter','latex','fontsize',24);
ylabel('I(V)','interpreter','latex','fontsize',24);
title(['I-V for $\frac{\Delta_1}{\Delta_2} = 1.0$'],'interpreter','latex','fontsize',20);

