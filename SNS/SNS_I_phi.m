% SNS Josephson Effect with the BdG Hamiltonain

% Calculation of I vs E 
% tight-binding parameter [eV]
t = 1.0;

% 0+ for the iteration to converge to one of the roots
eta = 1e-8;

% N_D : number of device points
N_D = 2;

% Electrochemical potentials of the contacts
% mu : Fermi level of the device
mu = 0.0;
mu1 = mu;
mu2 = mu;

% Temperature kT : [eV]
kT = 0.0001;

N_phi = 50;
phi_vec = 2*pi*linspace(0,1,N_phi);
I_phi_vec = zeros(1,length(phi_vec));

for ii = 1:length(phi_vec)
    
    % phase difference betweem the two superconductors
    phi = phi_vec(ii)
    
    % Superconducting order paramter
    Delta1 = 0.001;
    Delta2 = 0.001*exp(1j*phi);
    
    I_E = @(E) calculate_I_E(E,t,mu1,mu2,Delta1,Delta2,kT,eta,N_D);
    I_phi_vec(ii) = quadv(I_E,-4*Delta1,-Delta1,1e-4) + quadv(I_E,Delta1,4*Delta1,1e-4);
    
end


figure(1)
plot(phi_vec,I_phi_vec,'linewidth',2.0);
xlabel('$\phi$','interpreter','latex','fontsize',16);
ylabel('$I(\phi)$','interpreter','latex','fontsize',16)
title('Josephson Effect','fontsize',16,'interpreter','latex');