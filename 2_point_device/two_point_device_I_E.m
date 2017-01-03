% Modelling a 2-point channel connected to infinite contacts
% No BdG Hamiltonian, just simple TB model
% Author: Sandesh Kalantre

% Calculation of I vs E from the analytic result and comparison to purely
% numerical result

% tight-binding parameter [eV]
t = 1.0;

% 0+ for the iteration to converge to one of the roots
eta = 1e-8;

% Device Hamiltonian
H_D = [2*t -t;-t 2*t];

% Electrochemical potentials of the contacts
mu1 = 2+1;
mu2 = 2-1;

% Temperature kT : [eV]
kT = 0.001;

% Fermi function
fermi = @(E,mu,kT) 1.0/(1.0 + exp((E-mu)/kT));

% N_E : number of points in the energy grid
N_E = 100;
E_vec = linspace(-t,5*t,N_E); 

I_vec = zeros(1,length(E_vec));
I_vec_analytic = zeros(1,length(E_vec));

for ii = 1:length(E_vec)
    E = E_vec(ii);
    
    g1 = surface_g_numerical(E,t,0,eta);
    g2 = surface_g_numerical(E,t,0,eta);
    
    Sigma1 = [t*t*g1 0;0 0];
    Sigma2 = [0 0;0 t*t*g2];
    
    Gamma1 = 1j*(Sigma1 - Sigma1');
    Gamma2 = 1j*(Sigma2 - Sigma2');
    
    G_D = inv((E  +1j*eta)*eye(2) - H_D - Sigma1 - Sigma2);
    
    Sigma_corr = Gamma1*fermi(E,mu1,kT) + Gamma2*fermi(E,mu2,kT);
    G_corr = G_D * Sigma_corr * G_D';
    
    I_vec(ii) = 1j*(H_D(1,2)*G_corr(2,1) - G_corr(1,2)*H_D(2,1));
    
    %analytic expression for I(E)
    gamma1 = 1j*t^2*(g1 - g1');
    gamma2 = 1j*t^2*(g2 - g2');
    
    num = t^2*gamma1*gamma2;%*(fermi(E,mu2,kT) - fermi(E,mu1,kT));
    den = abs((E -2*t - t^2*g1)*(E -2*t - t^2*g2) - t^2)^2;
    I_vec_analytic(ii) = num/den;
    %I_vec_analytic(ii) = (fermi(E,mu2,kT) - fermi(E,mu1,kT));
end

figure(1)
%plot(E_vec,I_vec,'linewidth',2.0);
xlabel('$\frac{E}{t}$','interpreter','latex','fontsize',16);
ylabel('T(E)','fontsize',16)
title('Transmission Function','fontsize',16,'interpreter','latex');
hold on;
plot(E_vec,I_vec_analytic,'linewidth',2.0);
hold off;