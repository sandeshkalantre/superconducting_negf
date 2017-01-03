% DOS for a contact using the surface Green's function i.e local DOS at the
% interface

% Calculation of DOS vs E
clear;

% eta : 0+ for calculating retarded Green's functions
eta = 1e-10;
    
% t0 : tight binding parameter = hbar^2/(2 m a^2)
% t0 units : eV
t0 = 1.0;

% kT : k_B * Temperature(K)
% kT units : eV
% kT default value : 0.026 eV (Room temperature)
kT = 0.0001;

% Potential of the contact
mu = 1.0;

% Change in contact Hamiltonian because of bias
alpha = [2*t0 + mu  0;0 -2*t0 - mu];
beta = -t0* [1 0; 0 -1];

% N_E : number of points in the energy grid
N_E = 1000;
E_vec = 5* t0.* linspace(-1,1,N_E);

% local density of states
a = zeros(1,length(E_vec));

% electron and hole density
n = zeros(1,length(E_vec));


for ii = 1:length(E_vec)
    E = E_vec(ii);
    
    g = surface_g(E,alpha,beta,eta);
    A = 1j*(g - g');
    
    fermi = @(E,mu,kT) 1.0/(1.0 + exp((E - mu)/kT));
    Fermi = [fermi(E,-mu,kT) 0;0 1-fermi(E,mu,kT)];
    
    N = A * Fermi;
    a(ii) = trace(A);
    n(ii) = trace(N);
end

figure(1)
plot(E_vec,a,'linewidth',1.5);
xlabel('Energy');
ylabel('DOS');

figure(2)
plot(E_vec,n,'linewidth',1.5);
xlabel('Energy');
ylabel('e-/h density');
