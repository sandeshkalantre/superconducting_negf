function g = surface_g_analytic(E,t,mu)
    % Calculates the surface Green's function as a root of the quadratic
    % equation
    % Only one of the roots is returned, check the final code to see which
    % is actually returned
    
    % +/- sign before the square root 
    num = (E - mu -2*t) + sqrt((-E + mu +2*t)^2 - 4*t*t); 
    den = 2*t*t;
    g = num/den;
end