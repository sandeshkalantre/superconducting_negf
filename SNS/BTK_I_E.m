function I = BTK_I_E(E,mu1,mu2,Delta,kT)
fermi = @(E,mu,kT) 1.0/(1.0 + exp((E - mu)/kT));
V = mu1 - mu2;
if abs(E) < Delta
    I = 2 * (fermi(E,V,kT) - fermi(E,-V,kT));
else
    I = (1 + (Delta^2)/(abs(E) + sqrt(E^2 - Delta^2))^2) *  (fermi(E,V,kT) - fermi(E,-V,kT));
end

