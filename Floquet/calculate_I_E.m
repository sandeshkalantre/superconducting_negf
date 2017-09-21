function I_E = calcuate_I_E(E,t,trans,U,mu,mu1,mu2,Delta1,Delta2,kT,eta,N_D)
    % Physical constants
    q = 1.0;
    hbar = 1.0;
    

    % Device Hamiltonian (BdG Hamiltonian in the transformed domain)
    % We assume mu = 0 and Delta = 0 in the device region

    alpha = kron(eye(3),[2*t - mu  0; 0 -2*t + mu]);
    beta = kron(eye(3),-t* [1 0; 0 -1]);

    H_D = zeros(6*N_D,6*N_D);

    for jj = 1:N_D
        if jj == 1
            H_D(2*jj-1,2*jj-1) = alpha(1,1);
            H_D(2*jj-1,2*jj) = alpha(1,2) + Delta1;
            H_D(2*jj,2*jj-1) = alpha(2,1) + conj(Delta1);
            H_D(2*jj,2*jj) = alpha(2,2); 
        elseif jj == N_D
            H_D(2*jj-1,2*jj-1) = alpha(1,1);
            H_D(2*jj-1,2*jj) = alpha(1,2) + Delta2;
            H_D(2*jj,2*jj-1) = alpha(2,1) + conj(Delta2);
            H_D(2*jj,2*jj) = alpha(2,2);    
        else        
            H_D(2*jj-1,2*jj-1) = alpha(1,1);
            H_D(2*jj-1,2*jj) = alpha(1,2);
            H_D(2*jj,2*jj-1) = alpha(2,1);
            H_D(2*jj,2*jj) = alpha(2,2);
        end
        
        %off diagonal terms
        if(jj < N_D)
            H_D(2*jj-1,2*jj+1) = beta(1,1);
            H_D(2*jj-1,2*jj+2) = beta(1,2);
            H_D(2*jj,2*jj+1) = beta(2,1);
            H_D(2*jj,2*jj+2) = beta(2,2);

            H_D(2*jj+1,2*jj-1) = beta(1,1);
            H_D(2*jj+1,2*jj) = beta(1,2);
            H_D(2*jj+2,2*jj-1) = beta(2,1);
            H_D(2*jj+2,2*jj) = beta(2,2);
        end
    end

    alpha1 = kron(eye(3),[2*t - mu  Delta1; conj(Delta1) -2*t + mu ]);

    alpha2 = kron(eye(3),[2*t - mu   Delta2; conj(Delta2) -2*t + mu ]);
    
    g1 = surface_g(E,alpha1,beta,eta);
    g2 = surface_g(E,alpha2,beta,eta);
    
    Sigma1 = zeros(2*N_D);
    g1 = beta_trans*g1*beta_trans';
    Sigma1(1,1) = g1(1,1);
    Sigma1(1,2) = g1(1,2);
    Sigma1(2,1) = g1(2,1);
    Sigma1(2,2) = g1(2,2);
    Gamma1 = 1j*(Sigma1 - Sigma1');
    
    Sigma2 = zeros(2*N_D);
    g2 = beta_trans*g2*beta_trans';
    Sigma2(2*N_D - 1,2*N_D - 1) = g2(1,1);
    Sigma2(2*N_D - 1,2*N_D) = g2(1,2);
    Sigma2(2*N_D,2*N_D - 1) = g2(2,1);
    Sigma2(2*N_D,2*N_D) = g2(2,2);
    Gamma2 = 1j*(Sigma2 - Sigma2');
    
    G_D = inv((E + 1j*eta) .* eye(2*N_D) - H_D - Sigma1 - Sigma2);
    
    fermi = @(E,mu,kT) 1.0/(1.0 + exp((E - mu)/kT));
    
    Fermi_matrix1 = [fermi(E,mu1-mu,kT) 0;0 fermi(E,-mu1+mu,kT)];
    Fermi1 = kron(eye(N_D),Fermi_matrix1);
    
    Fermi_matrix2 = [fermi(E,mu2-mu,kT) 0;0 fermi(E,-mu2+mu,kT)];
    Fermi2 = kron(eye(N_D),Fermi_matrix2);
    
    Sigma_corr = Gamma1*Fermi1 + Gamma2*Fermi2;
    G_corr = G_D*Sigma_corr*G_D';
    
    %alternate way for G_corr
    %Sigma = Sigma1 + Sigma2;
    %g = inv((E + 1j*eta) .* eye(2*N_D) - H_D);
    %Fermi_matrix = [fermi(E + mu,mu1,kT) 0;0 fermi(-E - mu,-mu,kT)];
    %Fermi = kron(eye(N_D),Fermi_matrix);
    %a = 1j*(g - g') * Fermi;
    %G_corr = (eye(2*N_D) + Sigma*G_D)*a*(eye(2*N_D) + G_D'*Sigma');
    
    I_op = (1j*q/hbar)*(H_D(1:2,3:4)*G_corr(3:4,1:2) - H_D(3:4,1:2)*G_corr(1:2,3:4));
    I_E = 2*I_op(1,1); 
end