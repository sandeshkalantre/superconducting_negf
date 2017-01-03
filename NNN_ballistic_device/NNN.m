% NNN Ballistic Device
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
t0 = 1;

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

% N_V : number of points in the voltage grid
N_V = 50;
V_vec = 1 * linspace(-1,1,N_V);
I_vec = zeros(1,length(V_vec));

for ii = 1:length(V_vec)
    V = V_vec(ii)
    
    % mu1 : electrochemical potential of the left contact
    mu1 = 0.0 + V/2;
    % mu2 : electrochemical potential of the right contact
    mu2 = 0.0 - V/2;

    alpha1 = [2*t0 + mu1 0; 0 -2*t0 - mu1];

    alpha2 = [2*t0 + mu2 0; 0 -2*t0 - mu2];
    
    I = 0;
    
    % N_E : number of points in the energy grid for integration
    N_E = 1000;
    E_vec = 0.2 * t0 * linspace(-1,1,N_E);
    deltaE = E_vec(2) - E_vec(1);

    %integration using rectangle rule
    for jj = 1:length(E_vec)
        E = E_vec(jj);
        
        g1 = surface_g(E,alpha1,beta,eta);
        g2 = surface_g(E,alpha2,beta,eta);
        
        Sigma1 = zeros(2*N_D);
        g1 = beta*g1*beta';
        Sigma1(1,1) = g1(1,1);
        Sigma1(1,2) = g1(1,2);
        Sigma1(2,1) = g1(2,1);
        Sigma1(2,2) = g1(2,2);
        Gamma1 = 1j*(Sigma1 - Sigma1');
        
        Sigma2 = zeros(2*N_D);
        g2 = beta*g2*beta';
        Sigma2(2*N_D - 1,2*N_D - 1) = g2(1,1);
        Sigma2(2*N_D - 1,2*N_D) = g2(1,2);
        Sigma2(2*N_D,2*N_D - 1) = g2(2,1);
        Sigma2(2*N_D,2*N_D) = g2(2,2);
        Gamma2 = 1j*(Sigma2 - Sigma2');
        
        G_D = inv((E + 1j*eta) .* eye(2*N_D) - H_D - Sigma1 - Sigma2);
        
        fermi = @(E,mu,kT) 1.0/(1.0 + exp((E - mu)/kT));
        Fermi_matrix1 = [fermi(E,-mu1,kT) 0;0 1-fermi(E,mu1,kT)]; 
        Fermi1 = kron(eye(N_D),Fermi_matrix1);
        
        Fermi_matrix2 = [fermi(E,-mu2,kT) 0;0 1-fermi(E,mu2,kT)]; 
        Fermi2 = kron(eye(N_D),Fermi_matrix2);
        
        Sigma_corr = Gamma1*Fermi1 + Gamma1*Fermi2;
        G_corr = G_D*Sigma_corr*G_D';
        
        I_op = (1j*q/hbar) .* trace(G_D * Gamma1 * G_D' * Gamma2) * (fermi(E,mu1,kT) - fermi(E,mu2,kT));
        I = I + q*deltaE/(2*pi) * I_op;
    end
    I_vec(ii) = I;
end

figure(1)
plot(V_vec,I_vec,'linewidth',2.0);
set(gca,'FontSize',20)
xlabel('V','interpreter','latex','fontsize',20);
ylabel('Current (A)','interpreter','latex','fontsize',20);
title(['I-V for ' num2str(N_D) ' device point(s) at temperature ' num2str(kT) ' eV'],'interpreter','latex','fontsize',20);
