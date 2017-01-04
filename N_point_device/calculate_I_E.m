function I_E = calcuate_I_E(E,t,mu1,mu2,kT,eta,N_D)
    % Device Hamiltonian
    H_D = full(gallery('tridiag',N_D,-t,2*t,-t));

    % Fermi function
    fermi = @(E,mu,kT) 1.0/(1 + exp((E-mu)/kT));

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
    
    I_E = 1j*(H_D(1,2)*G_corr(2,1) -G_corr(1,2)*H_D(2,1));
end