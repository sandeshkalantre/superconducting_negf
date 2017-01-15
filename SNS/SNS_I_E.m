% Modelling a N-point channel connected to infinite contacts
% BdG Hamiltonian
% Author: Sandesh Kalantre

% Calculation of I vs E 
% tight-binding parameter [eV]
t = 1.0;

% 0+ for the iteration to converge to one of the roots
eta = 1e-6;

% N_D : number of device points
N_D = 100;

% Electrochemical potentials of the contacts
% mu : Fermi level of the device
mu = 0.0;
mu1 = mu;
mu2 = mu;

% Temperature kT : [eV]
kT = 0.0001;

% phase difference betweem the two superconductors
phi = 0;

% Superconducting order paramter
Delta1 = 0.001;
Delta2 = 0.001*exp(1j*phi);

% N_E : number of points in the energy grid
N_E = 100;
E_vec = 0.001* linspace(-1,1,N_E); 

I_vec = zeros(1,length(E_vec));

for ii = 1:length(E_vec)
    E = E_vec(ii);
    
    I_vec(ii) = calculate_I_E(E,t,mu1,mu2,Delta1,Delta2,kT,eta,N_D);
 end

figure(1)
plot(E_vec,I_vec,'linewidth',2.0);
xlabel('$\frac{E}{t}$','interpreter','latex','fontsize',16);
ylabel('I(E)','interpreter','latex','fontsize',16)
title('I(E)','fontsize',16,'interpreter','latex');
