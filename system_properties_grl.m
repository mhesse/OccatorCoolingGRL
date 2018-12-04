function [phi,kappa_sys,rhocp_sys,H,rho,cp,kappa,h,chi] = system_properties_grl(Ts,Tl,L,phi,sim)
% author: Marc A. Hesse
% date: 8 May 2018

% DESCRIPTION:
% Defines the physical properties of five phase system composing Ceres'
% crust.

% INPUT:
% Ts = solidus temperature [K]
% Tl = liquidus temperature [K]
% L = latent heat of water [J/kg]

% phi = structure containing the function handles for the volume fractions 
%       of the phases as function of temperature [-]
% rho = structure containing the constant densities of the phases [kg/m^3]
% cp = structure containing the function handles for the specific heat 
%      capacities of the phases as function of temperature [J/(kg K)]
% kappa = structure containing the function handles for the thermal 
%      conductivities of the phases as function of temperature [W/(m K)]
% h = structure containing the function handles for the specific enthalpy
%      of the phases as function of temperature [J/kg]
% chi = function handle of indicator function for melting interval [-]

if strcmp(sim.type,'grl')
    [phi,rho,cp,kappa,h,chi] = physical_properties_grl(Ts,Tl,L,phi);
elseif strcmp(sim.type,'linear_conduction')
    [phi,rho,cp,kappa,h,chi] = physical_properties_lin(Ts,Tl,L,phi);
elseif strcmp(sim.type,'Bowling')
    [phi,rho,cp,kappa,h,chi] = physical_properties_lin(Ts,Tl,L,phi);
else
end


DT = Tl - Ts;  % Melting interval

%% Total enthalpy of the system
H = @(T) rho.ice*phi.ice(T).*h.ice(T) + ...
         rho.hyd*phi.hyd(T).*h.hyd(T) + ...
         rho.sal*phi.sal(T).*h.sal(T) + ...
         rho.sil*phi.sil(T).*h.sil(T) + ...
         rho.wat*phi.wat(T).*h.wat(T);

kappa_sys = @(T) phi.ice(T).*kappa.ice(T) + ...
                 phi.hyd(T).*kappa.hyd(T) + ...
                 phi.sal(T).*kappa.sal(T) + ...
                 phi.sil(T).*kappa.sil(T) + ...
                 phi.wat(T).*kappa.wat(T);

%% Effective heat capacity of system
rhocp_mean = @(T) rho.ice*phi.ice(T).*cp.ice(T) + ...
                  rho.hyd*phi.hyd(T).*cp.hyd(T) + ...
                  rho.sal*phi.sal(T).*cp.sal(T) + ...
                  rho.sil*phi.sil(T).*cp.sil(T) + ...
                  rho.wat*phi.wat(T).*cp.wat(T);

latent_heat = @(T) chi(T)/DT.*(rho.wat*phi.wat_star*h.wat(T)...
                              -rho.ice*phi.ice0*h.ice(T)...
                              -rho.hyd*phi.hyd0*h.hyd(T)...
                              -rho.sal*phi.sal0*h.sal(T));             

rhocp_sys = @(T) rhocp_mean(T) + latent_heat(T);
end