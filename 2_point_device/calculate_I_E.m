function I_E = calcuate_I_E(E,t,mu1,mu2,kT,eta)
    % Device Hamiltonian
    H_D = [2*t -t;-t 2*t];
    
    % Fermi function
    fermi = @(E,mu,kT) 1.0/(1 + exp((E-mu)/kT));
           
    g1 = surface_g_numerical(E,t,0,eta);
    g2 = surface_g_numerical(E,t,0,eta);
    
    Sigma1 = [t*t*g1 0;0 0];
    Sigma2 = [0 0;0 t*t*g2];
    
    Gamma1 = 1j*(Sigma1 - Sigma1');
    Gamma2 = 1j*(Sigma2 - Sigma2');
    
    G_D = inv((E  +1j*eta)*eye(2) - H_D - Sigma1 - Sigma2);
    
    Sigma_corr = Gamma1*fermi(E,mu1,kT) + Gamma2*fermi(E,mu2,kT);
    G_corr = G_D * Sigma_corr * G_D';
    
    I_E = 1j*(H_D(1,2)*G_corr(2,1) -G_corr(1,2)*H_D(2,1));
    
    %analytic expression for I(E)
    %gamma1 = 1j*t^2*(g1 - g1');
    %gamma2 = 1j*t^2*(g2 - g2');
    
    %num = t^2*gamma1*gamma2*(fermi(E,-mu2,kT) - fermi(E,-mu1,kT));
    %den = abs((E -2*t - t^2*g1)*(E -2*t - t^2*g2) - t^2)^2;
    %I_E = num/den;
end