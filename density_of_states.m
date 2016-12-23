% Superconducting NEGF code for SNS junction
% Author: Sandesh Kalantre

% Current operator represented in the guage transformed domain
% Does not reproduce AC Josephson effects at finite bias

% DOS calculation including the effects of the leads
% Andreev bound states visible in the gap for N_D > 10

% eta : 0+ for calculating retarded Green's functions
eta = 1e-8;
    
% number of points in the device (channel)
N_D = 30;

% t0 : tight binding parameter = hbar^2/(2 m a^2)
% t0 units : eV
t0 = 1.0;

% Device Hamiltonian (BdG Hamiltonian in the transformed domain)
% We assume mu = 0 and Delta = 0 in the device region

alpha = [2*t0  0; 0 -2*t0];
beta = -t0* [1 0; 0 -1];

H_D = zeros(2*N_D,2*N_D);

for ii = 1:N_D
        H_D(2*ii-1,2*ii-1) = alpha(1,1);
        H_D(2*ii-1,2*ii) = alpha(1,2);
        H_D(2*ii,2*ii-1) = alpha(2,1);
        H_D(2*ii,2*ii) = alpha(2,2);

        %off diagonal terms
        if(ii < N_D)
            H_D(2*ii-1,2*ii+1) = beta(1,1);
            H_D(2*ii-1,2*ii+2) = beta(1,2);
            H_D(2*ii,2*ii+1) = beta(2,1);
            H_D(2*ii,2*ii+2) = beta(2,2);
       
            H_D(2*ii+1,2*ii-1) = beta(1,1);
            H_D(2*ii+1,2*ii) = beta(1,2);
            H_D(2*ii+2,2*ii-1) = beta(2,1);
            H_D(2*ii+2,2*ii) = beta(2,2);
       end    
end

% kT : k_B * Temperature(K)
% kT units : eV
% kT default value : 0.026 eV (Room temperature)
kT = 0.026;

% Delta: Superconducting order paramter
% Delta units : eV
Delta = 0.01;

% mu1 : electrochemical potential of the left contact
mu1 = 0.0;
% mu2 : electrochemical potential of the right contact
mu2 = 0.3;

% phi : phase difference between the two superconductors
phi = 0.5*pi;

Delta1 = Delta;
alpha1 = [2*t0 - mu1 Delta1; conj(Delta1) -2*t0 + mu1];
beta1 = -t0* [1 0; 0 -1];

Delta2 = Delta * exp(1j*phi);
alpha2 = [2*t0 - mu2 Delta2; conj(Delta2) -2*t0 + mu2];
beta2 = -t0* [1 0; 0 -1];

% N_E : number of points in the energy grid
N_E = 1000;
E_vec = 2 * Delta .* linspace(-1,1,N_E);

%density of states
A = zeros(1,length(E_vec));

for ii = 1:length(E_vec)
    E = E_vec(ii);

    g1 = surface_g(E,alpha1,beta1,eta);
    g2 = surface_g(E,alpha2,beta2,eta);
    
    Sigma1 = zeros(2*N_D);
    Sigma1(1,1) = g1(1,1);
    Sigma1(1,2) = g1(1,2);
    Sigma1(2,1) = g1(2,1);
    Sigma1(2,2) = g1(2,2);
    
    Sigma2 = zeros(2*N_D);
    Sigma2(2*N_D - 1,2*N_D - 1) = g1(1,1);
    Sigma2(2*N_D - 1,2*N_D) = g1(1,2);
    Sigma2(2*N_D,2*N_D - 1) = g1(2,1);
    Sigma2(2*N_D,2*N_D) = g1(2,2);
    
    G_D = inv((E + 1i*eta) .* eye(2*N_D) - H_D - Sigma1 - Sigma2);
    
    A(ii) = trace(1j * (G_D - G_D'));
end

figure(1)
plot(E_vec,A,'linewidth',1.5)
xlabel('Energy ($\frac{E}{\Delta}$)','interpreter','latex','fontsize',16);
ylabel('Density of States','interpreter','latex','fontsize',16);
title(['DOS for number of device points : ' num2str(N_D)],'interpreter','latex','fontsize',16)

