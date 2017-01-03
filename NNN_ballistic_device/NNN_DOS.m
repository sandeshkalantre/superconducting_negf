% NNN Ballistic Device
% Author: Sandesh Kalantre

% Calculation of DOS vs E
clear;
% kT : k_B * Temperature(K)
% kT units : eV
% kT default value : 0.026 eV (Room temperature)
kT = 0.0001;

[E_vec,A1,N1] = NNN_DOS_function(0,0,kT);
[E_vec,A2,N2] = NNN_DOS_function(0,1,kT);
[E_vec,A3,N3] = NNN_DOS_function(1,0,kT);

figure(1)
plot(E_vec,A1,'linewidth',2.0)
hold on;
plot(E_vec,A2,'linewidth',2.0)
plot(E_vec,A3,'linewidth',2.0)
hold off;

set(gca,'FontSize',20)
xlabel('Energy ($\frac{E}{t_0}$)','interpreter','latex','fontsize',20);
ylabel('Density of States','interpreter','latex','fontsize',20);
%ylim([100,101]);
legend('0 V' ,'1 V' ,'-1 V');

figure(2)
plot(E_vec,N1,'linewidth',1.5);
hold on;
plot(E_vec,N2,'linewidth',1.5);
plot(E_vec,N3,'linewidth',1.5);
hold off;
legend('0 V' ,'1 V' ,'-1 V');