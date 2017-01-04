% Modelling a N-point channel connected to infinite contacts
% No BdG Hamiltonian, just simple TB model
% Author: Sandesh Kalantre

% Calculation of I vs E from the analytic result
% tight-binding parameter [eV]
t = 1.0;

% 0+ for the iteration to converge to one of the roots
eta = 1e-8;

% Device Hamiltonian
% N_D : number of points in the device region
N_D = 2;
H_D = full(gallery('tridiag',N_D,-t,2*t,-t));

% Electrochemical potentials of the contacts
% mu : Fermi level of the device
mu = 2.0;
mu1 = mu + 0.5;
mu2 = mu - 0.5;

% Temperature kT : [eV]
kT = 0.001;

% Fermi function
fermi = @(E,mu,kT) 1.0/(1.0 + exp((E-mu)/kT));

% N_E : number of points in the energy grid
N_E = 100;
E_vec = linspace(-t,5*t,N_E); 

I_vec = zeros(1,length(E_vec));

for ii = 1:length(E_vec)
    E = E_vec(ii);
    
    g1 = surface_g_numerical(E,t,0,eta);
    g2 = surface_g_numerical(E,t,0,eta);
    
    Sigma1 = zeros(N_D);
    Sigma1(1,1) = t*t*g1;
    Sigma2 = zeros(N_D);
    Sigma2(N_D,N_D) = t*t*g2;
    
    Gamma1 = 1j*(Sigma1 - Sigma1');
    Gamma2 = 1j*(Sigma2 - Sigma2');
    
    G_D = inv((E  +1j*eta)*eye(N_D) - H_D - Sigma1 - Sigma2);
    
    Sigma_corr = Gamma1*fermi(E,mu1,kT) + Gamma2*fermi(E,mu2,kT);
    G_corr = G_D * Sigma_corr * G_D';
    
    I_vec(ii) = 1j*(H_D(1,2)*G_corr(2,1) - G_corr(1,2)*H_D(2,1));   
 end

figure(1)
plot(E_vec,I_vec,'linewidth',2.0);
xlabel('$\frac{E}{t}$','interpreter','latex','fontsize',16);
ylabel('I(E)','interpreter','latex','fontsize',16)
title('I(E)','fontsize',16,'interpreter','latex');
