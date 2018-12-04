function [] = plot_volume_fractions(Tmin,Tmax,Ts,Tl,phi,col)
T = [Tmin,Ts,Tl,Tmax];
Temp = [T,fliplr(T)];

%% Define Bottom and Top of areas
null = zeros(size(T));
sil_cum = phi.sil(T);
sal_cum = sil_cum+phi.sal(T);
hyd_cum = sal_cum+phi.hyd(T);
ice_cum = hyd_cum+phi.ice(T);
wat_cum = ice_cum+phi.wat(T);

patch(Temp,[null,fliplr(sil_cum)],col.tan,'LineStyle','none'), hold on
patch(Temp,[sil_cum,fliplr(sal_cum)],col.red,'LineStyle','none')
patch(Temp,[sal_cum,fliplr(hyd_cum)],col.green,'LineStyle','none')
patch(Temp,[hyd_cum,fliplr(ice_cum)],col.orange,'LineStyle','none')
patch(Temp,[ice_cum,fliplr(wat_cum)],col.blue,'LineStyle','none')

plot([Ts Ts],[0 1],'k--')
plot([Tl Tl],[0 1],'k--')
text(200,.83,'ICE','fontsize',14)
text(200,.55,'CLATHRATE','fontsize',14)
text(200,.1,'SILICATE ROCK','fontsize',14)
text(200,.3,'HYDRATED SALT','fontsize',14)
text(270,.7,'BRINE','fontsize',14)


end