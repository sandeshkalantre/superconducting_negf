% QP current using the tunnel Hamiltonian approach
% Author: Sandesh Kalantre

% The device does not have any points. The setup consists on two
% superconductors coupled by a coupling element M whose strength can be
% varied

% Current operator represented in the guage transformed domain
% Does not reproduce AC Josephson effects at finite bias

% Calculation of the current vs voltage relation

clear;

% eta : 0+ for calculating retarded Green's functions
eta = 1e-8;

% Fundamental Constants
q = 1.6e-19;
hbar = 1.1e-34;

% kT : k_B * Temperature(K)
% kT units : eV
% kT default value : 0.026 eV (Room temperature)
kT = 0.00001;

% t0 : tight binding parameter = hbar^2/(2 m a^2)
% t0 units : eV
t0 = 1.0;

% Delta: Superconducting order paramter
% Delta units : eV
Delta = 0.001;

% Coupling element for the Tunnel
M_strength = 0.001;
M = M_strength * [-1 0;0 1];

% phi : phase difference between the two superconductors
phi = 0;

% N_V : number of points in the voltage grid
N_V = 100;
V_vec = 3 * Delta .* linspace(-1,1,N_V);
I_vec = zeros(1,length(V_vec));

for ii = 1:length(V_vec)
    V = V_vec(ii)
    
    % mu1 : electrochemical potential of the left contact
    mu1 = 0.0 + V/2;
    % mu2 : electrochemical potential of the right contact
    mu2 = 0.0 - V/2;

    Delta1 = 0*Delta;
    alpha1 = [2*t0+mu1 Delta1; conj(Delta1) -2*t0-mu1];
    beta1 = -t0* [1 0; 0 -1];

    Delta2 = Delta * exp(1j*phi);
    alpha2 = [2*t0+mu2 Delta2; conj(Delta2) -2*t0-mu2];
    beta2 = -t0* [1 0; 0 -1];
    
    I = 0;
    
    % N_E : number of points in the energy grid for integration
    N_E = 200;
    E_vec = 10 * Delta .* linspace(-1,1,N_E);
    deltaE = E_vec(2) - E_vec(1);
    
    %integration using rectangle rule
    A1 = zeros(1,length(E_vec));
    A2 = zeros(1,length(E_vec));
    fermi_func = zeros(1,length(E_vec));
    for jj = 1:length(E_vec)
        E = E_vec(jj);
        
        g1 = surface_g_tunnel(E,alpha1,beta1,eta);
        g2 = surface_g_tunnel(E,alpha2,beta2,eta);
        
        %G = inv([inv(g1) -M;-M' inv(g2)]);
        fermi = @(E,mu,kT) 1.0/(1.0 + exp((E - mu)/kT));
        
        %sigma1 = beta1*g1*beta1';
        %sigma2 = beta2*g2*beta2';
        %gamma1 = 1j.*(sigma1 - sigma1');
        %gamma2 = 1j.*(sigma2 - sigma2');
        %Sigma_corr = [gamma1*fermi(E,mu1,kT) zeros(2);zeros(2) gamma2*fermi(E,mu2,kT)];
        %G_corr = G*Sigma_corr*G';
        
        %A =  1j.*(G - G');
        a1 = 1j*(g1 - g1');
        a2 = 1j*(g2 - g2');
        
        A1(jj) = trace(a1);
        A2(jj) = trace(a2);
        %Gamma1 = [gamma1 zeros(2);zeros(2) zeros(2)];
        %Gamma2 = [zeros(2) zeros(2);zeros(2) gamma2];
        %I_op = (q/hbar) * (G*Gamma1*G'*Gamma2) * (fermi(E,mu1,kT) - fermi(E,mu2,kT));
        %I_op = (1j*q/hbar).*(M * G_corr(3:4,1:2) - G_corr(1:2,3:4)*M');
        
    
        I = I + q*deltaE * trace(a1)*trace(a2) * (fermi(E,-mu1,kT) + fermi(E,-mu2,kT));
        fermi_func(jj) = fermi(E,mu1,kT) - fermi(E+mu2,0,kT);
    end
    %figure(2)
    %plot(E_vec,A1);
    %hold on;
    %plot(E_vec,A2);
    %plot(E_vec,fermi_func);
    %hold off;
    
    I_vec(ii) = I;
end

clf;
figure(1)
plot(V_vec./Delta,I_vec,'linewidth',2.0)
%hold on;
%plot(V_vec./Delta,imag(I_vec),'linewidth',2.0)
%hold off;
xlabel('$\frac{eV}{\Delta}$','interpreter','latex','fontsize',24);
ylabel('Current (arb units)','interpreter','latex','fontsize',24);
title(['Current-Voltage relation for a tunnel Josephson Junction'] ,'interpreter','latex','fontsize',24);
str = {['$t_0 =  $ ' num2str(t0) ' eV '],['$\Delta = $ ' num2str(Delta) ' eV'], ['$M = $ ' num2str(M_strength) ' eV']};
dim = [0.6 0.6 0.2 0.2];
annotation('textbox',dim,'String',str,'FitBoxToText','on','interpreter','latex','fontsize',24);