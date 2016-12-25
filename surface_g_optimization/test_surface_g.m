% Optimization of the chi value (see definition of surface_g) to improve
% speed of the fixed point iteration
% Author: Sandesh Kalantre

eta = 1e-8;
tolerance = 1e-2;
t0 = 1.0;
Delta = 0.01;
alpha = [2*t0 Delta;conj(Delta) -2*t0];
beta = -t0 * [1 0;0 -1];

N_chi = 100;
chi_vec = linspace(0,2,N_chi);
iter_count_vec = zeros(1,length(chi_vec));
E = 0.0;

for ii = 1:length(chi_vec)
    chi = chi_vec(ii)
    [g , iter_count] = surface_g(E,alpha,beta,eta,chi,tolerance); 
    iter_count_vec(ii) = iter_count;
end

figure(1)
semilogy(chi_vec(2:length(chi_vec)),iter_count_vec(2:length(iter_count_vec)),'linewidth',1.5);
xlabel('$\chi$','interpreter','latex');
ylabel('Number of iterations')
title(['Energy ' num2str(E) ' | Tolerance: ' num2str(tolerance)])