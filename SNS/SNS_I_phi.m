% SNS Josephson Effect with the BdG Hamiltonain
clear;
% Calculation of I vs E 
% tight-binding parameter [eV]
t = 1.0;

% transmission paramter
% t value for the coupling element
trans = 1.0

% 0+ for the iteration to converge to one of the roots
eta = 1e-8;

% N_D : number of device points
N_D = 50;

% Electrochemical potentials of the contacts
% mu : Fermi level of the device
mu = 2.0;
mu1 = mu;
mu2 = mu;

% Temperature kT : [eV]
kT = 0.0001;

N_phi = 50;
phi_vec = 2*pi*linspace(0,1,N_phi);
I_phi_vec = zeros(2,length(phi_vec));

for ii = 1:length(phi_vec)
    
    % phase difference betweem the two superconductors
    phi = phi_vec(ii)
    
    U = 0;
    
    % Superconducting order paramter
    Delta1 = 0.01;
    Delta2 = 0.01*exp(1j*phi);
    I_E = @(E) calculate_I_E(E,t,trans,U,mu,mu1,mu2,Delta1,Delta2,kT,eta,N_D); 
    I_phi_vec(1,ii) = integral(I_E,-1.5*Delta1,1.5*Delta1,'AbsTol',1e-8,'ArrayValued',true);
    
    Delta1 = 0.005;
    Delta2 = 0.005*exp(1j*phi);
    I_E = @(E) calculate_I_E(E,t,trans,U,mu,mu1,mu2,Delta1,Delta2,kT,eta,N_D);
    I_phi_vec(2,ii) = integral(I_E,-1.5*Delta1,1.5*Delta1,'AbsTol',1e-8,'ArrayValued',true);

end

figure(1)
plot(phi_vec,I_phi_vec(1,:),'linewidth',2.5);
hold on;
plot(phi_vec,I_phi_vec(2,:),'linewidth',2.5);
hold off;
xlabel('$\phi$','interpreter','latex','fontsize',24);
ylabel('$I(\phi)$','interpreter','latex','fontsize',24)
title('Josephson Effect','fontsize',24,'interpreter','latex');




