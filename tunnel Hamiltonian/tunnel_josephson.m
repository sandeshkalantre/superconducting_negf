% Josephson Effect using the tunnel Hamiltonian approach
% Author: Sandesh Kalantre

% The device does not have any points. The setup consists on two
% superconductors coupled by a coupling element M whose strength can be
% varied

% Current operator represented in the guage transformed domain
% Does not reproduce AC Josephson effects at finite bias

% Calculation of the current vs phase relation

clear;

% eta : 0+ for calculating retarded Green's functions
eta = 1e-8;

% Fundamental Constants
q = 1.6e-19;
hbar = 1.1e-34;

% kT : k_B * Temperature(K)
% kT units : eV
% kT default value : 0.026 eV (Room temperature)
kT = 0.0001;

% t0 : tight binding parameter = hbar^2/(2 m a^2)
% t0 units : eV
t0 = 1.0;

% Delta: Superconducting order paramter
% Delta units : eV
Delta = 0.001;

% mu1 : electrochemical potential of the left contact
mu1 = 0.0;
% mu2 : electrochemical potential of the right contact
mu2 = 0.0;

% Coupling element for the Tunnel
M_strength = 0.01;
M = M_strength * [-1 0;0 1];

% N_phi : number of points in the phi grid
N_phi = 50;
phi_vec = 2*pi*linspace(0,1,N_phi);

I_vec = zeros(1,length(phi_vec));

for ii = 1:length(phi_vec)
    phi = phi_vec(ii)
    
    Delta1 = Delta;
    alpha1 = [2*t0 - mu1 Delta1; conj(Delta1) -2*t0 + mu1];
    beta1 = -t0* [1 0; 0 -1];

    Delta2 = Delta * exp(1j*phi);
    alpha2 = [2*t0 - mu2 Delta2; conj(Delta2) -2*t0 + mu2];
    beta2 = -t0* [1 0; 0 -1];
    
    I = 0;
    
    % N_E : number of points in the energy grid for integration
    N_E = 25;
    E_vec = 10 * Delta .* linspace(-1,1,N_E);
    deltaE = E_vec(2) - E_vec(1);
    
    for jj = 1:length(E_vec)
        E = E_vec(jj);
        
        g1 = surface_g(E,alpha1,beta1,eta);
        g2 = surface_g(E,alpha2,beta2,eta);
        
        G = inv([inv(g1) -M;-M' inv(g2)]);
        fermi = @(E,mu,kT) 1.0/(1.0 + exp((E - mu)/kT));
        
        sigma1 = beta1*g1*beta1';
        sigma2 = beta2*g2*beta2';
        gamma1 = 1j*(sigma1 - sigma1');
        gamma2 = 1j*(sigma2 - sigma2');
        Sigma_corr = [gamma1*fermi(E,mu1,kT) zeros(2);zeros(2) gamma2*fermi(E,mu2,kT)];
        G_corr = G*Sigma_corr*G';
        I_op = (1j*q/hbar) * (G_corr(3:4,1:2)*M' -  M * G_corr(1:2,3:4));
        
        I = I + q*deltaE * (I_op(1,1) - I_op(2,2));
    end
    
    I_vec(ii) = I;
end

clf;
figure(1)
plot(phi_vec./(2*pi),I_vec,'linewidth',2.0)
xlabel('Phase $\frac{\phi}{2 \pi}$','interpreter','latex','fontsize',24);
ylabel('Current (arb units)','interpreter','latex','fontsize',24);
title(['Current-phase relation for a tunnel Josephson Junction'] ,'interpreter','latex','fontsize',24);
str = {['$t_0 =  $ ' num2str(t0) ' eV '],['$\Delta = $ ' num2str(Delta) ' eV'], ['$M = $ ' num2str(M_strength) ' eV']};
dim = [0.6 0.6 0.2 0.2];
annotation('textbox',dim,'String',str,'FitBoxToText','on','interpreter','latex','fontsize',24);