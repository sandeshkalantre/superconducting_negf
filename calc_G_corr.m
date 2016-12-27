function G_corr = calc_G_corr(E,H_D,g1,g2,t0,mu1,mu2,kT,eta)
    N_D = 0.5 * length(H_D);
    
    beta = -t0 * [1 0;0 -1];

    Sigma1 = zeros(2*N_D);
    g1 = beta * g1 *beta';
    Sigma1(1,1) = g1(1,1);
    Sigma1(1,2) = g1(1,2);
    Sigma1(2,1) = g1(2,1);
    Sigma1(2,2) = g1(2,2);

    Gamma1 = 1j*(Sigma1 - Sigma1');

    Sigma2 = zeros(2*N_D);
    g2 = beta * g2 *beta';
    Sigma2(2*N_D - 1,2*N_D - 1) = g2(1,1);
    Sigma2(2*N_D - 1,2*N_D) = g2(1,2);
    Sigma2(2*N_D,2*N_D - 1) = g2(2,1);
    Sigma2(2*N_D,2*N_D) = g2(2,2);

    Gamma2 = 1j*(Sigma2 - Sigma2');

    G_D = inv((E + 1i*eta) .* eye(2*N_D) - H_D - Sigma1 - Sigma2);
    fermi = @(E,mu,kT) 1.0/(1.0 + exp((E - mu)/kT));
    
    Sigma_corr = Gamma1*fermi(E,mu1,kT) + Gamma2*fermi(E,mu2,kT);
    G_corr = G_D*Sigma_corr*G_D';
end