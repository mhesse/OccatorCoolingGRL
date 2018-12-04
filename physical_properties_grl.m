function [phi,rho,cp,kappa,h,chi] = physical_properties_grl(Ts,Tl,L,phi)
% author: Marc A. Hesse
% date: 7 May 2018

% DESCRIPTION:
% Defines the physical properties of the phases in Ceres' volatile rich
% crust.

% INPUT:
% Ts = solidus temperature [K]
% Tl = liquidus temperature [K]
% L = latent heat of water [J/kg]
% phi = structure containing the initial volume fractions of phases [-]
%       phi.ice0, phi.hyd0, phi.sal0, phi.sil0

% OUTPUT:
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

% EXAMPLE CALL
% >> Ts = 245;
% >> Tl = 273;
% >> L = 3.34e5;
% >> phi0.ice = 0.3;
% >> phi0.hyd = 0.3;
% >> phi0.sal = 0.2;
% >> phi0.sil = 0.2;
% >> [phi,rho,cp,kappa,h,chi] = physical_properties_grl(Ts,Tl,L,phi0);

%% Indicator function for the melting inteval [-]
% This is simply a function that is 1 for Ts < T < Tl and 0 elsewhere
chi = @(T) (T>Ts).*(T<Tl);

%% Densities of the phases [kg/m^3]
rho.ice = 917;
rho.hyd = 1000;
rho.sal = 2200;
rho.sil = 2430;
rho.wat = 1000;

%% Volume fractions of the phases as function of T [-]
f = @(T) (Tl-T)/(Tl-Ts); % function decreasing linearly from 1 to 0 between Ts and Tl
phi.wat_star = phi.ice0+phi.hyd0+phi.sal0;
phi.ice = @(T) phi.ice0*(T<=Ts) + phi.ice0*f(T).*chi(T); 
phi.hyd = @(T) phi.hyd0*(T<=Ts) + phi.hyd0*f(T).*chi(T); 
phi.sal = @(T) phi.sal0*(T<=Ts) + phi.sal0*f(T).*chi(T); 
phi.sil = @(T) phi.sil0 + T*0;
phi.wat = @(T) phi.wat_star*(T>=Tl)+phi.wat_star*(1-f(T)).*chi(T);

%% Heat capacity [J/kg/K]
% All phases considered here have at most a linear variation in cp
heat_capacity = @(T,const) const.a + const.b*T;

% Coefficients in the cp polynomial
cp.const.ice.a = 185;  cp.const.ice.b = 7.037;
cp.const.hyd.a = 494;  cp.const.hyd.b = 6.1;
cp.const.wat.a = 4200; cp.const.wat.b = 0;
cp.const.sil.a = 2000; cp.const.sil.b = 0;
cp.const.sal.a = 920;  cp.const.sal.b = 0;

cp.ice = @(T) heat_capacity(T,cp.const.ice);
cp.hyd = @(T) heat_capacity(T,cp.const.hyd);
cp.wat = @(T) heat_capacity(T,cp.const.wat);
cp.sil = @(T) heat_capacity(T,cp.const.sil);
cp.sal = @(T) heat_capacity(T,cp.const.sal);

%% Thermal conductivity [W/m/K]
kappa.ice = @(T) 0.4685 + 488.12./T;
kappa.hyd = @(T) 0.64 + T*0;
kappa.sal = @(T) 0.6  + T*0;
kappa.sil = @(T) 2 + T*0;
kappa.wat = @(T) 0.56 + T*0;

%% Specific enthalpy of the phases
% General integral of the linear heat capacity function
integral_heat_capacity = @(T,const) const.a*(T-Ts) + const.b/2 * (T.^2 - Ts^2);
h.ice = @(T) integral_heat_capacity(T,cp.const.ice);
h.hyd = @(T) integral_heat_capacity(T,cp.const.hyd);
h.sal = @(T) integral_heat_capacity(T,cp.const.sal);
h.sil = @(T) integral_heat_capacity(T,cp.const.sil);
h.wat = @(T) integral_heat_capacity(T,cp.const.wat) + L;

