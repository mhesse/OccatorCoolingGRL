function [uss,Lss,rss] = compute_geotherm(Grid,Kd,Tbot,Ttop,D,G,I,fs,B,N,fn,Param)

% Residual of steady-state heat equation
Lss = @(u) -D*Kd(u)*G;
rss = @(u) Lss(u)*u-fs;

% Initial guess is linear interpolation between Ttop and Tbot
u = solve_lbvp(-D*G,0+fn,B,Param.g,N);

nres = norm(rss(u)); 
ndu = 1; i = 0; tol = 1e-6; imax = 15;
fprintf('Solving for steady-state geotherm:\n')
while (nres > tol || ndu > tol) && i < imax
    uold = u;
    u = solve_lbvp(Lss(u),fs+fn,B,Param.g,N);
    res = rss(u); res(Param.dof_dir) = 0; res(Param.dof_neu) = 0;
    nres = norm(res(:)); ndu = norm(u-uold); i = i+1;
    fprintf('it = %d: nres = %3.2e ndu = %3.2e.\n',i,nres,ndu)
    
end

if i >= imax
    error('Steady-state geotherm did not converge in %d iterations!',imax)
else
    fprintf('Steady-state geotherm converged to tol = %3.2e.\n',tol)
    uss = u;
end
