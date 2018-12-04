function [umax,zmax,umin,zmin,vol] = solution_post_processing(u,uss,R,Z,Grid,Ts,phi)
% Note: The simple indexing below works because of y-first numbering!
%% Find maximum temperature along centerline
% Subtracting geotherm to find maximum in the interior.
du_center = u(Grid.dof_xmin) - uss(Grid.dof_xmin);

[dumax,imax] = max(du_center);
umax = u(imax); zmax = Z(imax); 

% this typically puts the max above real max
while imax > 1 && u(imax-1) > umax 
    umax = u(imax-1); imax = imax-1;
end

if imax > 1
    %% Find the minimum temperature below umax along center line
    ulower = u(1:imax);
    [umin,imin] = min(ulower); zmin = Z(imin);
    
    %% Find proper maximum
    u_center =  u(Grid.dof_xmin); z_center = Z(Grid.dof_xmin);
    [umax,imax] = max(u_center(imin:end)); zmax = z_center(imin-1+imax);
    
    %% Determine volume
    dof_chamber = Grid.dof( Z(:)>=zmin & phi.wat(u)>0);
    dof_partial = dof_chamber(phi.wat(u(dof_chamber)) < 1-phi.sil0);
    vol.total = phi.wat(u(dof_chamber))'*Grid.V(dof_chamber);
    vol.partial = phi.wat(u(dof_partial))'*Grid.V(dof_partial);
else
    fprintf('No more local maximum in T.\nStopping simulation.\n')
    umax = nan; umin = nan; 
    zmax = nan; zmin = nan;
    vol.total = nan; vol.partial = nan;
end

end