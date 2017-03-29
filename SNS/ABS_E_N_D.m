%spectrum 2D plot as a function of phase and N_D


% Superconducting order paramter
Delta1 = 1e-2;
Delta2 = 1e-2*exp(1j*phi);

% eta : 0+ for calculating retarded Green's functions
eta = 1e-8;

% Fundamental Constants
q = 1;
hbar = 1;

% t0 : tight binding parameter = hbar^2/(2 m a^2)
% t0 units : eV
t = 1;

% transmission paramter
% t value for the coupling element
trans = 1;

% mu : Device Fermi Level
mu = 2;

% kT : k_B * Temperature(K)
% kT units : eV
% kT default value : 0.026 eV (Room temperature)
kT = 0.00001;

phi = pi/2;
N_E = 1000;
N_D_max = 51;
%DOS = zeros(N_D,N_E);

E_vec = Delta1*linspace(-1,1,N_E);
N_D_vec = 2:1:N_D_max;
DOS = [];
for ii = 1:length(N_D_vec)
    N_D = N_D_vec(ii);
    DOS_E = @(E) calculate_DOS(E,t,trans,U,mu,mu1,mu2,Delta1,Delta2,kT,eta,N_D);
    DOS = [DOS ;arrayfun(DOS_E,E_vec)];
end

figure(1)
surf(E_vec,N_D_vec,DOS,'EdgeColor','None');

    
        