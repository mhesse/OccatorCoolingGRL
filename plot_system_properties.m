function [] = plot_system_properties(phi,Ts,Tl,rhocp,kappa,col)
Tmin = 150;
Tmax = 300;
tol = 1e-4;
Tplot = [linspace(Tmin,Ts-tol,30),...
       linspace(Ts+tol,Tl-tol,30),...
       linspace(Tl+tol,Tmax,30)];

figure('name','Overview of system properties')
subplot 311
plot_volume_fractions(Tmin,Tmax,Ts,Tl,phi,col)
ylabel('$\phi_p$','interp','latex','fontsize',16)
set(gca,'xtick',[150,200,Ts,250,Tl,300],'xticklabel',{'','','T_s','','T_l',''})
set(gca,'ytick',[0:.2:1])
set(gca,'fontsize',14)

subplot 312
plot(Tplot,rhocp(Tplot)/1e6,'k-','linewidth',1.5), hold on
plot([Ts Ts],[0 16],'k--')
plot([Tl Tl],[0 16],'k--')
ylim([0 16])
set(gca,'xtick',[150,200,Ts,250,Tl,300],'xticklabel',{'','','T_s','','T_l',''})
set(gca,'ytick',[0:4:16])
set(gca,'fontsize',14)
ylabel('$\widetilde{\rho c_p}$ [MJ/(m$^3$K)]','interp','latex','fontsize',16)

subplot 313
plot(Tplot,kappa(Tplot),'k-','linewidth',1.5), hold on
plot([Ts Ts],[0 3],'k--')
plot([Tl Tl],[0 3],'k--')
ylim([0 3])
set(gca,'xtick',[150,200,Ts,250,Tl,300],'xticklabel',{'150','200','','250','','300'})
set(gca,'fontsize',14)
xlabel('temperature [K]','interp','latex','fontsize',16)
ylabel('$\overline{\kappa}$ [W/m/K]','interp','latex','fontsize',16)
drawnow
end