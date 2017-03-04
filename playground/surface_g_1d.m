% Plot of real and imaginary parts of surface Green's function
clear;

% tight-binding parameter
t = 1;
% Fermi level
mu = 2;

N_E = 100;
E_vec = 4*t*linspace(-1,1,N_E);
g = zeros(2,length(E_vec));

for ii = 1:length(E_vec)
    E = E_vec(ii);
    g(1,ii) = ((E - 2*t + mu) + sqrt((E - 2*t + mu)^2 - 4*t^2))/(2*t^2);
    g(2,ii) = ((E - 2*t + mu) - sqrt((E - 2*t + mu)^2 - 4*t^2))/(2*t^2);    
end

figure(1)
plot(E_vec,real(g(1,:)));
hold on;
plot(E_vec,imag(g(1,:)));
hold off;
legend('Real','Imaginary');