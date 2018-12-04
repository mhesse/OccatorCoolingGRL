function [T] = initial_condition_grl(Tini,Tss,Geom,R,Z,Grid,Tl,Ts)
T = Tss;

%% Cryomagma chamber
rc = Geom.rc;
rh = Geom.rh;
zc = Grid.ymax-Geom.dc;
zh = Grid.ymax-Geom.dh;
T(R(:) <= rc & Z(:) >=zc) = Tini;

%% Mush halo in horizontal direction
T_hor = Tl - (Tl-Ts)/(Geom.rh-Geom.rc) * (R - Geom.rc); 
dof_hor = Grid.dof(R(:) > rc & R(:) <= rh & Z(:) >=zc);
T(dof_hor) = T_hor(dof_hor);

%% Mush halo in vertical direction
T_ver = Ts + (Tl-Ts)/(zc-zh) * (Z - zh); 
dof_ver = Grid.dof(Z(:) > zh & Z(:) <= zc & R(:) <= rc);
T(dof_ver) = T_ver(dof_ver);

%% Mush halo corner
A = [1 rc zc;...
     1 rc zh;...
     1 rh zc];
 b = [Tl;Ts;Ts];
x = A\b;
T_cor = x(1) + x(2)*R + x(3)*Z;

dof_cor = Grid.dof(R(:)<rh & R(:) > rc & Z(:) > zh & Z(:) < zc & T_cor(:) >= Ts);

T(dof_cor) = T_cor(dof_cor);