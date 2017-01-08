% Modelling a N-point channel connected to infinite contacts
% BdG Hamiltonian
% Author: Sandesh Kalantre

% Calculation of I vs E from the analytic result
% tight-binding parameter [eV]
t = 1.0;

% 0+ for the iteration to converge to one of the roots
eta = 1e-8;

% N_D : number of device points
N_D = 2;

% Electrochemical potentials of the contacts
% mu : Fermi level of the device
mu = 2.0;
mu1 = mu - 0.05;
mu2 = mu + 0.05;

% Temperature kT : [eV]
kT = 0.001;


% N_E : number of points in the energy grid
N_E = 100;
E_vec = linspace(-5*t,5*t,N_E); 

I_vec = zeros(1,length(E_vec));

for ii = 1:length(E_vec)
    E = E_vec(ii);
    
    I_vec(ii) = calculate_I_E(E,t,mu1,mu2,kT,eta,N_D);
 end

figure(1)
plot(E_vec,I_vec,'linewidth',2.0);
xlabel('$\frac{E}{t}$','interpreter','latex','fontsize',16);
ylabel('I(E)','interpreter','latex','fontsize',16)
title('I(E)','fontsize',16,'interpreter','latex');
