function DOS = calculate_DOS(E,t,trans,U,mu,mu1,mu2,Delta1,Delta2,kT,eta,N_D)
    % Physical constants
    q = 1.0;
    hbar = 1.0;
    
    
    % Device Hamiltonian (BdG Hamiltonian in the transformed domain)
    % We assume mu = 0 and Delta = 0 in the device region

    alpha = [2*t - mu  0; 0 -2*t + mu];
    beta = -t* [1 0; 0 -1];
    beta_trans = -trans* [1 0; 0 -1];
    
    H_D = zeros(2*N_D,2*N_D);

    for jj = 1:N_D
        if jj == 1
            H_D(2*jj-1,2*jj-1) = alpha(1,1) + U;
            H_D(2*jj-1,2*jj) = alpha(1,2) + Delta1;
            H_D(2*jj,2*jj-1) = alpha(2,1) + conj(Delta1);
            H_D(2*jj,2*jj) = alpha(2,2) - U; 
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

    alpha1 = [2*t - mu  Delta1; conj(Delta1) -2*t + mu ];

    alpha2 = [2*t - mu   Delta2; conj(Delta2) -2*t + mu ];
    
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
    A = 1j*(G_D - G_D');
    DOS = trace(A);
    