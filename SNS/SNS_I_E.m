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
mu = 2;
mu1 = mu - 1e-4;
mu2 = mu;

% Temperature kT : [eV]
kT = 0.000001;

% phase difference betweem the two superconductors
phi = 0;

% Superconducting order paramter
Delta1 = 0e-3;
Delta2 = 1e-5*exp(1j*phi);

% N_E : number of points in the energy grid
N_E = 100;
%E_vec = linspace(-5e-3,-2e-3,100);
E_vec =  2e-3*linspace(0,1,N_E); 

I_vec = zeros(1,length(E_vec));
BTK_I = zeros(1,length(E_vec));

for ii = 1:length(E_vec)
    E = E_vec(ii)
    
    I_vec(ii) = calculate_I_E(E,t,mu,mu1,mu2,Delta1,Delta2,kT,eta,N_D);
    %BTK_I(ii) = BTK_I_E(E,mu1,mu2,Delta2,kT);
end

figure(1)
plot(E_vec,I_vec,'linewidth',2.0);
hold on;
%plot(E_vec,BTK_I,'linewidth',2.0);
hold off;
xlabel('$E$','interpreter','latex','fontsize',16);
ylabel('I(E)','interpreter','latex','fontsize',16)
title('I(E)','fontsize',16,'interpreter','latex');
