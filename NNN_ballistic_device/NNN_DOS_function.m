function [E_vec,A,N] = NNN_DOS_function(mu1,mu2,kT)
    % eta : 0+ for calculating retarded Green's functions
    eta = 1e-8;
    
    % number of points in the device (channel)
    N_D = 2;
    
    % t0 : tight binding parameter = hbar^2/(2 m a^2)
    % t0 units : eV
    t0 = 1.0;
    
    % Device Hamiltonian
    alpha = [2*t0  0; 0 -2*t0];
    beta = -t0* [1 0; 0 -1];
    
    H_D = zeros(2*N_D,2*N_D);
    
    for jj = 1:N_D
        H_D(2*jj-1,2*jj-1) = alpha(1,1);
        H_D(2*jj-1,2*jj) = alpha(1,2);
        H_D(2*jj,2*jj-1) = alpha(2,1);
        H_D(2*jj,2*jj) = alpha(2,2);
        
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

    % Change in contact Hamiltonian because of bias
    alpha1 = [2*t0 + mu1  0;0 -2*t0 - mu1];
    alpha2 = [2*t0 + mu2  0;0 -2*t0 - mu2];
    
    % N_E : number of points in the energy grid
    N_E = 1000;
    E_vec = 5* t0.* linspace(-1,1,N_E);
    
    % density of states
    A = zeros(1,length(E_vec));
    %e/h density
    N = zeros(1,length(E_vec));
    
    for ii = 1:length(E_vec)
        E = E_vec(ii);
    
        g1 = surface_g(E,alpha1,beta,eta);
        g2 = surface_g(E,alpha2,beta,eta);
        
        Sigma1 = zeros(2*N_D);
        Sigma1(1,1) = g1(1,1);
        Sigma1(1,2) = g1(1,2);
        Sigma1(2,1) = g1(2,1);
        Sigma1(2,2) = g1(2,2);
        Gamma1 = 1j*(Sigma1 - Sigma1');
        
        Sigma2 = zeros(2*N_D);
        Sigma2(2*N_D - 1,2*N_D - 1) = g2(1,1);
        Sigma2(2*N_D - 1,2*N_D) = g2(1,2);
        Sigma2(2*N_D,2*N_D - 1) = g2(2,1);
        Sigma2(2*N_D,2*N_D) = g2(2,2);
        Gamma2 = 1j*(Sigma2 - Sigma2');
        
        G_D = inv((E + 1j*eta) .* eye(2*N_D) - H_D - Sigma1 - Sigma2);
        
        fermi = @(E,mu,kT) 1.0/(1.0 + exp((E - mu)/kT));
        
        Fermi_matrix1 = [fermi(E,-mu1,kT) 0;0 1-fermi(E,mu1,kT)]; 
        Fermi1 = kron(eye(N_D),Fermi_matrix1);
        
        Fermi_matrix2 = [fermi(E,-mu2,kT) 0;0 1-fermi(E,mu2,kT)]; 
        Fermi2 = kron(eye(N_D),Fermi_matrix2);
        
        Sigma_corr = Gamma1*Fermi1 + Gamma2*Fermi2;
        G_corr = G_D*Sigma_corr*G_D';
        
        A_matrix = 1j*(G_D - G_D');
        
        A(ii) = trace(A_matrix);
        N(ii) = trace(G_corr);
    end
end