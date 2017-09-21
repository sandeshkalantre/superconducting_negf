% Point contact S1-N-S2 device
% The contacts are assumed to be a single point. The N region between the
% contacts has N_D device points.

% The entire device is characterised by a single tight-binding parameter t

% Each contact is characterised by an electrochemical potential mu and
% superconducting order paramter Delta

% [Physical Parameters]
% 0+ for the iteration to converge to the retarted Green's functions
eta = 1e-3;

% t : tight-binding parameter
t = 1.0;

% kT : k_B * Temperature(K)
% kT units : eV
% kT default value : 0.026 eV (Room temperature)
kT = 0.00001;

% [Device Parameters]
% N_D : number of device points
N_D = 2;

% mu : Device Fermi Level
mu = 2;

% [Contact Paramters]
% phi : phase difference betweem the two superconductors
phi = 0;

% superconducting order paramters for the contacts
Delta1 = 1e-3;
Delta2 = 0e-2*exp(1j*phi);

% electrochemical potential of the contacts
mu1 = mu + 1e-3;
mu2 = mu;

% N_E : number of points in the energy grid
N_E = 100;
E_vec =  2e-3*linspace(-1,1,N_E); 

I_vec = zeros(1,length(E_vec));

for ii = 1:length(E_vec)
    E = E_vec(ii)
    
    I_vec(ii) = calculate_I_E(E,t,mu,mu1,mu2,Delta1,Delta2,kT,eta,N_D);
end

figure(1)
plot(E_vec,I_vec,'linewidth',2.0);
xlabel('$E$','interpreter','latex','fontsize',16);
ylabel('I(E)','interpreter','latex','fontsize',16)
title('I(E)','fontsize',16,'interpreter','latex');


