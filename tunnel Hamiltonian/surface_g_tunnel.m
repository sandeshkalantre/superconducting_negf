function g = surface_g_tunnel(E,alpha,beta,eta)
    % error tolerance in the result
    tolerance = 1e-4;
    % maximum number of iterations
    N_lim = 50000;
 
    % error
    % set to 1 in at the start of the first iteration
    err = 1;
    
    % iter_count : keep track of number of iterations
    % if iter_counts exceeds N_lim, display 'Convergence failed' and value
    % of E
    iter_count = 0;
    
    % value of g set to inv(alpha) at the start of first iteration
    % it is possible to set other default values to achieve better
    % convergence
    g = inv(alpha);
    g_last = inv(alpha);
    
    while err > tolerance
        g = inv((E + 1j*eta).*eye(2) - alpha - beta'*g*beta);
        err = norm(g - g_last,1)/norm(g,1);
        
        if  err < tolerance
            break;
        end
        
        % for faster convergence
        % change 0.5 to different values for getting better results
        g = g_last + 0.5 * (g - g_last);
        g_last = g;
    
        iter_count = iter_count + 1;
        if(iter_count == N_lim)
            disp('Convergence failed for [E err]')
            disp([E err])
            break;
        end
    end 
end