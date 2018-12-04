function [results] = OccatorCoolingGRL(Tini,phi_hyd,chamber_size,col,fig_name,bottom_bc,time_save,sim)
% file: OccatorGRL_dimensional.m
% author: Marc A. Hesse
% date: 11 May 2018

% Description:
% Models the cooling of the presumed impact induced cryomagma chamber beneath
% Occator crater on Ceres. This code is dimensional!

% INPUT:
% Tini = Initial temperature of the fully molten region in the core of the 
%        chamber [K].
% phi_hyd = Initial hydrate volume fraction in the crust [-]
% chamber_size = sacling factor for the initial temperature distribution
% [-]

% OUTPUT:
% t_final = Time when the simulation was terminated because Tmax < Ts
%           this corresponds to the maximum age of the magma chamber.
% dH_total = ???
% Tlow = maximum temperature in the lower crust beneath the chamber.

%% Unit conversions
s2yr = 60^2*24*365.25;
results.s2Ma = s2yr*1e6;

%% Physical properties
xmax = 60e3;       % [m] Maximum radius of 92 km
ymin = -1e3;      % [m] Depth into the mantle
ymax = 46e3;       % [m] Depth of the crust
tmax = 20e6*s2yr;  % [s] Simulation time

Geom.rc = chamber_size*4e3;        % [m] Initial radius of chamber
Geom.dc = chamber_size*2.e4;       % [m] Initial depth of chamber
Geom.rh = chamber_size*1.2e4;      % [m] Initial radius of mushy halo
Geom.dh = chamber_size*2.5e4;      % [m] Initial depth of mushy halo
results.Geom = Geom;

Tl = 273.15;     % [K] liquidus
Ts = 245;        % [K] solidus
Ttop = 150;      % [K] surface temperature beneath thermal skin
Tbot = 250;      % [K] temperature at the base of the crust
if strcmp(sim.type,'grl')
    L = 3.34e5;      % [J/kg] Latent heat of water
elseif strcmp(sim.type,'linear_conduction') || strcmp(sim.type,'Bowling')
    L = 0;
else
    error('Undefined simualtion type.')
end

% Initial volume fractions
phi.sil0 = 0.20;
phi.sal0 = 0.15;
phi.hyd0 = phi_hyd;
phi.ice0 = 1-phi.sil0-phi.sal0-phi_hyd; % volume fraction constraint

% Create function handles for the thermo-physical properties
[phi,kappa_sys,rhocp_sys,H,rho,cp,kappa,h,chi] = system_properties_grl(Ts,Tl,L,phi,sim);
plot_system_properties(phi,Ts,Tl,rhocp_sys,kappa_sys,col)

%% Build Grid and Discrete Operators
%% Numerical Properties
Nx = round(xmax/250);
Ny = round((ymax-ymin)/250);
Nt = 650;

% Build Grid
Grid.xmin = 0; Grid.xmax = xmax; Grid.Nx = Nx; % vertical direction
Grid.ymin = ymin; Grid.ymax = ymax; Grid.Ny = Ny; % radial direction
if strcmp(sim.type,'grl') || strcmp(sim.type,'Bowling')
    Grid.geom = 'cylindrical_rz'; 
elseif strcmp(sim.type,'linear_conduction'); Grid.geom = 'cartesian'; end
Grid = build_grid(Grid);
[R,Z] = meshgrid(Grid.xc,Grid.yc);
% Grid.dof_int = setdiff(Grid.dof,Grid.dof_ymin);

% Build Operators
[D,G,I] = build_ops(Grid);
Kd = @(u) comp_mean(reshape(kappa_sys(u),Grid.Ny,Grid.Nx),-1,1,Grid); 
fs = sparse(Grid.N,1,0);

% Set boundary conditions


%% Solve for steady-state geotherm and build initial condition
% Due to non-linearity we first solve the Dirichlet problem
% and compute the appropriate basal heat flux from it. Then
% we set this as a Neumann boundary condition.
if strcmp(sim.type,'grl')
    % 1) Compute steady-state state
    if bottom_bc == 'neu'
        % Set boundary conditions
        Param.dof_dir   = Grid.dof_ymax;
        Param.dof_f_dir = Grid.dof_f_ymax;
        Param.g = Ttop*ones(Grid.Nx,1);
        
        Param.dof_neu   = Grid.dof_ymin;
        Param.dof_f_neu = Grid.dof_f_ymin;
        Param.qb = comp_basal_heat_flow(phi,Grid); % look-up precomputed heat flow
        
        % Solve steady state
        [B,N,fn] = build_bnd(Param,Grid,I);
        [uss,Lss,rss] = compute_geotherm(Grid,Kd,Tbot,Ttop,D,G,I,fs,B,N,fn,Param);
        
        % Extract basal temperature
        T_bot = interp1(Grid.yc,uss(Grid.dof_xmin),0);
        fprintf('Basal Temperature is %3.2f [K]\n',T_bot)
        results.T_bot = T_bot;
    else
        % Set boundary conditions
        Param.dof_dir   = [Grid.dof_ymin;  Grid.dof_ymax];
        Param.dof_f_dir = [Grid.dof_f_ymin;Grid.dof_f_ymax];
        Param.g = [Tbot*ones(Grid.Nx,1);Ttop*ones(Grid.Nx,1)];
        
        Param.dof_neu = [];
        Param.dof_f_neu = [];
        Param.qb = [];
        
        % Solve steady state
        [B,N,fn] = build_bnd(Param,Grid,I);
        [uss,Lss,rss] = compute_geotherm(Grid,Kd,Tbot,Ttop,D,G,I,fs,B,N,fn,Param);
        
        % Compute basal heat flow
        qss = comp_flux(D,Kd(uss),G,uss,fs,Grid,Param);
        fprintf('Basal heat flux is %3.2f [mW/m^2]\n',mean(qss(Grid.dof_f_ymin)))
        
    end
elseif strcmp(sim.type,'linear_conduction')
    % Set boundary conditions - no-flow everywhere!
    Param.dof_dir   = [];
    Param.dof_f_dir = [];
    Param.g = [];
    
    Param.dof_neu = [];
    Param.dof_f_neu = [];
    Param.qb = [];
    
    % Solve steady state
    [B,N,fn] = build_bnd(Param,Grid,I);
    Lss = @(u) -D*Kd(u)*G;
    rss = @(u) Lss(u)*u-fs;
elseif strcmp(sim.type,'Bowling')
    % Set boundary conditions
    Param.dof_dir   = [Grid.dof_ymin;  Grid.dof_ymax];
    Param.dof_f_dir = [Grid.dof_f_ymin;Grid.dof_f_ymax];
    Tbot = Ttop+0.5*Grid.Ly/1e3
    Param.g = [Tbot*ones(Grid.Nx,1);Ttop*ones(Grid.Nx,1)];
    
    Param.dof_neu = [];
    Param.dof_f_neu = [];
    Param.qb = [];
    
    % Solve steady state
    [B,N,fn] = build_bnd(Param,Grid,I);
    [uss,Lss,rss] = compute_geotherm(Grid,Kd,Tbot,Ttop,D,G,I,fs,B,N,fn,Param);
    
    % Compute basal heat flow
    qss = comp_flux(D,Kd(uss),G,uss,fs,Grid,Param);
    fprintf('Basal heat flux is %3.2f [mW/m^2]\n',mean(qss(Grid.dof_f_ymin)))
else
    error('Unknown simulation type.\n')
end

% figure
% surf(R/1e3,Z/1e3,reshape(uss,Grid.Ny,Grid.Nx))
% shading interp
% view(2)
% colorbar
% [uss,Lss,rss] = compute_geotherm(Grid,Kd,Tbot,Ttop,D,G,I,fs,B,N,fn,Param);
% norm(uss-uss_dir)

% 2) Add impact-induced cryomagma chamber
if strcmp(sim.type,'grl') || strcmp(sim.type,'Bowling')
    u0 = initial_condition_grl(Tini,uss,Geom,R,Z,Grid,Tl,Ts);
elseif strcmp(sim.type,'linear_conduction')
    r0 = 5e3; % initial radius
    uss = zeros(Grid.N,1); u0 = uss;
    u0(R(:)<=r0) = Tini;
else
    error('Unknown simulation type.\n')
end
energy = plot_initial_condition_grl(u0,uss,phi,H,Grid,R,Z);
results.energy = energy;
results.Ts = Ts; results.Tl = Tl;

% 3) Make folder and save initial condition
if length(time_save) > 0
    fprintf('\nSaving initial condition to file: results_0\n')
    
    save('results_0','u0','uss','R','Z','phi','kappa','cp','rho','energy','h','H','Grid','Tl','Ts','Tini','Geom');
end
% fprintf('Saving initial condition:'), tic;
% folder_name = ['dc',num2str(Geom.dc/1e3),'_dh',num2str(Geom.dh/1e3),'_hyd',num2str(phi.hyd0*100),'_Tini',num2str(Tini)];
% cd ../../../Data/GRL/Dimensional
% mkdir(folder_name), cd(folder_name)
% save('time_0','u0','uss','R','Z','phi','kappa','cp','rho','energy','h','H','Grid','Tl','Ts','Tini','Geom');
% fprintf(' %3.2e s.\n',toc)

%% Solve transient problem
% Transient operators and residual
M = @(u) spdiags(rhocp_sys(u),0,Grid.N,Grid.N);
Im = @(u,dt) M(u) + dt*Lss(u);
Ex = @(u,dt) M(u);
r =  @(u,uold,dt) M(u)*(u-uold) + dt*rss(u);

tvec = union(tmax*(linspace(0,1,Nt)'.^2),time_save*results.s2Ma); 
Nt = length(tvec);

dtvec = [0;diff(tvec)];
umax = zeros(Nt,1); 
umean = umax;
vol.total = umax;
vol.partial = umax;
umin = umax;
[umax(1),zmax,umin(1),zmin,vol_temp] = solution_post_processing(u0,uss,R,Z,Grid,Ts,phi);
vol.total(1) = vol_temp.total; results.V = vol_temp.total;
vol.partial(1) = vol_temp.partial;
u = u0; time  = 0;
tol = 1e-4; imax = 6;
total_solves = 0;
results.connection = 0;
figure('name',fig_name)
t = 2;
while t <= Nt
    uold = u; dt = dtvec(t); time = tvec(t);
    fprintf('\n\n%d) time = %3.3f ka dt = %3.2e ka.\n',t-1,time/s2yr/1e3,dt/s2yr/1e3);
    nres = 1; ndu = 1; i = 0; nres_temp = 2;
    while (nres > tol || ndu > tol) && (i < imax && nres <= nres_temp)
        utemp = u; i = i+1; nres_temp = nres;
        u = solve_lbvp(Im(u,dt),Ex(u,dt)*uold+dt*(fs+fn),B,Param.g,N); 
        total_solves = total_solves+1;
        res = r(u,uold,dt); 
        res(Param.dof_dir) = 0; 
        res(Param.dof_neu) = 0; 
        nres = norm(res(:)); 
        du = u - utemp; ndu = norm(du); 
        if i == 1; 
            nres0 = nres; ndu0 = ndu;
            fprintf('it = 0: nres = %3.2e ndu = %3.2e.\n',nres0,ndu0)
        end
        if nres0 > 1e1; nres = nres/nres0; end;
        if ndu0 > 1e1;  ndu  = ndu/ndu0; end;
        fprintf('it = %d: nres = %3.2e ndu = %3.2e.\n',i,nres,ndu)
    end
    u = print_convergence_message(i,imax,nres,nres0,nres_temp,tol,u,utemp);
    q = comp_flux(D,Kd(u),G,u,fs,Grid,Param);
    [umax(t),zmax,umin(t),zmin,vol_temp] = solution_post_processing(u,uss,R,Z,Grid,Ts,phi);
    vol.total(t) = vol_temp.total;
    vol.partial(t) = vol_temp.partial;
    if strcmp(sim.type,'grl') || strcmp(sim.type,'Bowling')
        plot_transient_solution(R,Z,u,umax,umin,zmax,zmin,tvec,tmax,s2yr,Grid,t,Ts,Tl,vol,q,kappa_sys(u))
    elseif strcmp(sim.type,'linear_conduction')
        subplot 121
        Dsys = kappa_sys(250)/rhocp_sys(250);
        eta = R(Grid.dof_ymin)/sqrt(4*Dsys*tvec(t));
        theta = rhocp_sys(250)/abs(energy.H0_total)*sqrt(4*Dsys*tvec(t))*u(Grid.dof_ymin);
        plot(eta,theta,'r')
        subplot 122
        loglog(tvec(t)/s2yr/1e6,u(Grid.dof_xmin),'ro'), hold on
        drawnow
    end
    [save_true,save_index] = ismember(time,time_save*results.s2Ma);
    
    if save_true
        file_save = ['results_',num2str(save_index)];
        fprintf('\nSaving solution at time %3.2f [Ma] to file: %s\n',time/results.s2Ma,file_save)
        u_max = umax(t); u_min = umin(t); vol_tot = vol_temp.total; vol_par = vol_temp.partial;
        save(file_save,'u','time','u_max','u_min','vol_tot','vol_par')
    end
    
    if umin(t) > Ts 
        fprintf('Connection to deep reservoir Tmin = %3.2f K!\n',umin(t))
        results.connection = 1;
    end
    if ~strcmp(sim.type,'linear_conduction') && (isnan(umax(t)) || (umax(t) - umin(t) < 1e-3 && umax(t) < Ts))
        Nt = t;  % terminate the simulation
    end
    t = t+1;
end
fprintf('\nSimulation finished after %d timesteps (total linear solves = %d)\n',Nt,total_solves)
results.num.Nt = Nt;
results.num.lin_sol = total_solves;

results = comp_characteristic_quantities(results,tvec,umax,umin,vol,Nt);

% cd ../../../../MatlabCode/GRL/DimensionalCode

end

function [u] = print_convergence_message(i,imax,nres,nres0,nres_temp,tol,u,utemp)
if i >= imax+1 && nres <= nres_temp
        fprintf('Iteration did NOT converge in %d iterations!\n',imax)
        if nres0 > 1e1;
            fprintf('Final residual %3.2f %% of inital.\n',nres*100)
        else
            fprintf('Final residual %3.2f %% of inital.\n',nres/nres0*100)
        end
    elseif i < imax && nres > nres_temp
        fprintf('Residual stopped decreasing!\n')
        fprintf('Last iteration discarded.\n')
        if nres0 > 1e1;
            fprintf('Final residual %3.2f %% of inital.\n',nres*100)
        else
            fprintf('Final residual %3.2f %% of inital.\n',nres/nres0*100)
        end
        % Reset solution to previous iterate
        u = utemp;
    else
        fprintf('Iteration converged to tol = %3.2e.\n',tol)
    end
end



function [results] = comp_characteristic_quantities(results,tvec,umax,umin,vol,Nt)
%% Post-processing of the solution
tvec = tvec(umax ~= 0); 
umin = umin(umax~=0);
vol.total = vol.total(umax~=0);
vol.partial = vol.partial(umax~=0);
umax = umax(umax~=0); 

if isnan(umax(results.num.Nt))
    vol.total = vol.total(~isnan(umax));
    vol.partial = vol.partial(~isnan(umax));
    tvec = tvec(~isnan(umax)); 
    umin = umin(~isnan(umax)); 
    umax = umax(~isnan(umax)); 
end 
results.vol = vol;
results.tvec = tvec;
results.umin = umin;
results.umax = umax;
results.t_final = max(tvec);

% options = optimset('Display','iter'); % show iterations
options = optimset('Display','off');
% 1) Determine the cooling time max(T) = Ts
if umax(results.num.Nt-1) > results.Ts
    fprintf('Distinct pluton disappears before it cools below Ts.\n')
    results.t_cooling = tvec(Nt-1);
    fprintf('Final simulation time, %3.2f Ma, is set at cooling time.\n',results.t_cooling/results.s2Ma)
else
    obj_umax = @(t_int) abs(results.Ts - interp1(tvec,umax,t_int,'pchip'));
    t0_min = tvec(umax > results.Ts); t0_min = max(t0_min); % lower bound time
    t0_max = tvec(umax < results.Ts); t0_max = min(t0_max); % upper bound time
    t_cooling = fminbnd(obj_umax,t0_min,t0_max,options);
    results.t_cooling = t_cooling;
    fprintf('Pluton cools below Ts at %3.2f Ma.\n',t_cooling/results.s2Ma)
end


% 2) Determine the max T in the lower crust    

[umin0,imin] = min(-umin);
obj_umin = @(t_int) interp1(tvec,-umin,t_int,'pchip');
t_umin_max = fminbnd(obj_umin,tvec(imin-2),tvec(imin+2),options);
results.t_umin_max = t_umin_max;
results.umin_max = -obj_umin(results.t_umin_max);

% 3) Determine time of connection if relevant
if results.umin_max >= results.Ts
    fprintf('Connection to deep reservir occured\n')
    umin_con = umin(tvec<=results.t_umin_max);
    tvec_con = tvec(tvec<=results.t_umin_max);
    obj_con = @(t_int) abs(results.Ts - interp1(tvec_con,umin_con,t_int,'pchip'));
    t_connection = fminbnd(obj_con,0,t_umin_max,options);
    results.t_connection = t_connection;
end

end
