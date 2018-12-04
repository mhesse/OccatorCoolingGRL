function [] = plot_transient_solution(R,Z,u,umax,umin,zmax,zmin,tvec,tmax,s2yr,Grid,t,Ts,Tl,vol,q,kappa)
subplot 231
surf(R/1e3,Z/1e3,reshape(u,Grid.Ny,Grid.Nx)), view(2), shading interp
axis equal tight

subplot 232
plot(u(Grid.dof_xmin),Grid.yc(Grid.dof_xmin)/1e3), hold on
if umax(t) > Ts
    plot(umax(t),zmax/1e3,'ro')
else
    plot(umax(t),zmax/1e3,'ro','markerfacecolor','r')
end
if umin(t) < Ts
    plot(umin(t),zmin/1e3,'bo')
else
    plot(umin(t),zmin/1e3,'bo','markerfacecolor','b')
end
plot([Ts Ts],[Grid.ymin Grid.ymax]/1e3,'k--')
plot([Tl Tl],[Grid.ymin Grid.ymax]/1e3,'k--'), hold off
pbaspect([.3 1 1])

subplot 233
plot(tvec(1:t)/s2yr/1e6,umax(1:t)), hold on
plot(tvec(1:t)/s2yr/1e6,umin(1:t))
plot([0 tmax]/s2yr/1e6,[Ts Ts],'k--')
plot([0 tmax]/s2yr/1e6,[Tl Tl],'k--'), hold off
xlim([0,tmax]/s2yr/1e6), ylim([150 umax(1)+20])
xlabel 'time', ylabel 'T [K]'
pbaspect([1 .8 1])

subplot 234
plot(tvec(1:t)/s2yr/1e6,vol.total(1:t)/1e9), hold on
plot(tvec(1:t)/s2yr/1e6,vol.partial(1:t)/1e9), hold off
xlim([0,tmax]/s2yr/1e6), ylim([0 vol.total(1)/1e9])
xlabel 'time', ylabel 'V [km^3]'
pbaspect([1 .8 1])
drawnow


subplot 235
plot(Grid.xc/1e3,q(Grid.dof_f_ymin)*1e3,'-')%, hold on
% plot(Grid.xc/1e3,q(Grid.dof_f_ymax)*1e3,'--'), hold off
xlabel('Radius [km]')
ylabel('Heat flow [mW/m^2]')
ylim([0 8]), xlim([0 60])

subplot 236
plot(Grid.xc/1e3,kappa(Grid.dof_ymin),'-')%, hold on

xlabel('Radius [km]')
ylabel('Conductivity [W/m/K]')
ylim([0 4]), xlim([0 60])