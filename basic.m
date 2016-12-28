% Superconducting NEGF code for SNS junction
% Author: Sandesh Kalantre

% Current operator represented in the guage transformed domain
% Does not reproduce AC Josephson effects at finite bias

% Calculation of the current

% eta : 0+ for calculating retarded Green's functions
eta = 1e-8;

% Fundamental Constants
q = 1.6e-19;
hbar = 1.1e-34;
    
% number of points in the device (channel)
N_D = 1;

% t0 : tight binding parameter = hbar^2/(2 m a^2)
% t0 units : eV
t0 = 1.0;

% Device Hamiltonian (BdG Hamiltonian in the transformed domain)
% We assume mu = 0 and Delta = 0 in the device region

alpha = [2*t0  0; 0 -2*t0];
beta = -t0* [1 0; 0 -1];

H_D = zeros(2*N_D,2*N_D);

for jj = 1:N_D
        H_D(2*jj-1,2*jj-1) = alpha(1,1);
        H_D(2*jj-1,2*jj) = alpha(1,2);
        H_D(2*jj,2*jj-1) = alpha(2,1);
        H_D(2*jj,2*jj) = alpha(2,2);

        %off diagonal terms
        if(jj < N_D)
            H_D(2*jj-1,2*jj+1) = beta(1,1);
            H_D(2*jj-1,2*jj+2) = beta(1,2);
            H_D(2*jj,2*jj+1) = beta(2,1);
            H_D(2*jj,2*jj+2) = beta(2,2);
       
            H_D(2*jj+1,2*jj-1) = beta(1,1);
            H_D(2*jj+1,2*jj) = beta(1,2);
            H_D(2*jj+2,2*jj-1) = beta(2,1);
            H_D(2*jj+2,2*jj) = beta(2,2);
       end    
end

% kT : k_B * Temperature(K)
% kT units : eV
% kT default value : 0.026 eV (Room temperature)
kT = 0.0001;

% Delta: Superconducting order paramter
% Delta units : eV
Delta = 0.01;

% phi : phase difference between the two superconductors
phi = 0;

% N_V : number of points in the voltage grid
N_V = 50;
V_vec = 5 * Delta .* linspace(-1,1,N_V);
I_vec = zeros(1,length(V_vec));

for ii = 1:length(V_vec)
    V = V_vec(ii)
    
    % mu1 : electrochemical potential of the left contact
    mu1 = 0.0 +  V/2;
    % mu2 : electrochemical potential of the right contact
    mu2 = 0.0 - V/2;

    Delta1 = Delta;
    alpha1 = [2*t0 + mu1 Delta1; conj(Delta1) -2*t0 - mu1];
    beta1 = -t0* [1 0; 0 -1];

    Delta2 = Delta * exp(1j*phi);
    alpha2 = [2*t0 + mu2 Delta2; conj(Delta2) -2*t0 - mu2];
    beta2 = -t0* [1 0; 0 -1];
    
    I = 0;
    
    % N_E : number of points in the energy grid for integration
    N_E = 500;
    E_vec = 10 * Delta .* linspace(-1,1,N_E);
    deltaE = E_vec(2) - E_vec(1);
    
    %integration using rectangle rule
    for jj = 1:length(E_vec)
        E = E_vec(jj);
        
        g1 = surface_g(E,alpha1,beta1,eta);
        g2 = surface_g(E,alpha2,beta2,eta);
        
        fermi = @(E,mu,kT) 1.0/(1.0 + exp((E - mu)/kT));
        I = I + deltaE * (q*q/hbar) * trans(E,H_D,g1,g2,t0,eta) * (fermi(E,mu1,kT) - fermi(E,mu2,kT));
    end
    
    I_vec(ii) = I;
end

figure(1)
plot(V_vec./Delta,I_vec,'linewidth',2.0);
set(gca,'FontSize',20)
xlabel('$\frac{eV}{\Delta}$','interpreter','latex','fontsize',24);
ylabel('Current (A)','interpreter','latex','fontsize',24);
title(['I-V for ' num2str(N_D) ' device point(s) at temperature ' num2str(kT) ' eV'],'interpreter','latex','fontsize',24);
