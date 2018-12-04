function [q_bot] = comp_basal_heat_flow(phi,Grid)
load basal_heat_flow.mat

if Grid.Ny ~= Ny_basal + 4*abs(Grid.ymin)/1e3; error('Grid size does not match!\\ Re-run basal heatflow with your grid spacing'); end

q_bot = interp1(hyd_vec,qss_vec,phi.hyd0,'linear');  
    