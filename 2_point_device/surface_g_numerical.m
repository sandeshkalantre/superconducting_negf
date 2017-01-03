function g = surface_g_numerical(E,t,mu,eta)
    % Calculates the surgace Green's function using a fixed point iteration
    % scheme
    
    % final tolerance in the result
    tolerance = 1e-5;
    
    g = inv(E + 1j*eta - 2*t - mu);
    g_last = inv(E + 1j*eta - 2*t - mu);
    err = 1;
    while err > tolerance
        g = inv(E + 1j*eta - 2*t - mu - t*t*g);
        
        err = abs(g-g_last)/abs(g_last);
        
        g = g_last + 0.5*(g - g_last);
        g_last = g;
    end
end

