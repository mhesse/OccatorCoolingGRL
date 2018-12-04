function [energy] = plot_initial_condition_grl(u0,uss,phi,H,Grid,R,Z)

% Plot initial condition
figure('name','Initial Condition')
subplot 131
surf(R/1e3,Z/1e3,reshape(u0,Grid.Ny,Grid.Nx)), shading interp, view(2)
hold on
xlabel('radius [km]'), ylabel('z [km]')
title('Initial Temperature')
colorbar
axis equal
xlim([0 Grid.xmax/1e3]), ylim([0 Grid.ymax/1e3])

subplot 132
surf(R/1e3,Z/1e3,phi.wat(reshape(u0,Grid.Ny,Grid.Nx))), shading interp, view(2)
xlabel('radius [km]'), ylabel('z [km]')
title('Initial Melt Fraction')
colorbar
axis equal
xlim([0 Grid.xmax/1e3]), ylim([0 Grid.ymax/1e3])

subplot 133
H0 = H(reshape(u0,Grid.Ny,Grid.Nx));
Hss = H(reshape(uss,Grid.Ny,Grid.Nx));
surf(R/1e3,Z/1e3,(H0-Hss)/1e6), shading interp, view(2)
xlabel('radius [km]'), ylabel('z [km]')
title('Energy Added [MJ/m^3]')
colorbar
axis equal
xlim([0 Grid.xmax/1e3]), ylim([0 Grid.ymax/1e3])


energy.Hss_total = H(uss)'*Grid.V;
energy.H0_total = H(u0)'*Grid.V;
energy.dH_total = energy.H0_total - energy.Hss_total;

