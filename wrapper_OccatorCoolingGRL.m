%% Benchmark the OccatorGRL code vs linear Gaussian decay
% date: 11/12/2018
% author: Marc Hesse
close all, clear all, clc
% settings
sim.type = 'grl';
col = marc_colors;
Tini = 320;
phi_hyd = 0.35; % this will not matter because all properties are set to be same
chamber_size = 1.0;
bottom_bc = 'neu';

time_save = [.1:.1:4]; %[.1:.1:1.5]; % [Ma]
fig_name = 'Transient';
results = OccatorCoolingGRL(Tini,phi_hyd,chamber_size,col,fig_name,bottom_bc,time_save,sim);

save('Bowling_grl.mat','results','Tini','chamber_size','sim')