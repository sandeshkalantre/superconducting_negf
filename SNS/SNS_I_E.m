% Modelling a N-point channel connected to infinite contacts
% BdG Hamiltonian
% 
% Calculation of I vs E 
% tight-binding parameter [eV]
t = 1.0;

% 0+ for the iteration to converge to one of the roots
eta = 1e-8;

% N_D : number of device points
N_D = 2;

% Electrochemical potentials of the contacts
% mu : Fermi level of the device
mu = 0;
mu1 = mu;
mu2 = mu + 0.5e-2;

% Temperature kT : [eV]
kT = 0.000001;

% phase difference betweem the two superconductors
phi = 0;

% Superconducting order paramter
Delta1 = 1e-3;
Delta2 = 0*exp(1j*phi);

% N_E : number of points in the energy grid
N_E = 100;
%E_vec = linspace(-5e-3,-2e-3,100);
E_vec = 1e-2* linspace(-1,1,N_E); 
dE = E_vec(2)-E_vec(1);
I = 0;

I_vec = zeros(1,length(E_vec));

for ii = 1:length(E_vec)
    E = E_vec(ii);
    
    I_vec(ii) = calculate_I_E(E,t,mu1,mu2,Delta1,Delta2,kT,eta,N_D);
    I = I + I_vec(ii)*dE;
end
I

figure(1)
plot(E_vec,I_vec,'linewidth',2.0);
xlabel('$\frac{E}{t}$','interpreter','latex','fontsize',16);
ylabel('I(E)','interpreter','latex','fontsize',16)
title('I(E)','fontsize',16,'interpreter','latex');
