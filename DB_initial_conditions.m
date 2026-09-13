function [r, p] = DB_initial_conditions(obj)
% ==========================================================================
% Determines Computation of initial conditions depending on input
% - spherical: ICs for spherical orbits at arbitraty inclination
% - post-spherical: ICs for plunging quasi-spherical orbits (include small initial radial mom)
% - eccentric: ICs for eccentric orbits, both geodesic and non. Still preliminary
% Always assumes phi = 0 and pth = 0
% ==========================================================================
id = obj.ICs;

r = [obj.r0.*sin(obj.th0).*cos(obj.phi0);0;obj.r0.*cos(obj.th0)]; %3D orbits

if strcmp(id, 'spherical')
    % spherical initial conditions
    py0 = DB_spherical_ICs(obj);
    if obj.DBvKOS~=0
        py0 = obj.DBvKOS./obj.r0;
    end
    %py0 = 3.3742717620867./obj.r0; %specific values for elliptic orbits (used to compare with KOS)
    %py0 = 4.03253243659069./obj.r0;
    p   = [0;py0;0];

elseif strcmp(id, 'post-spherical')
    p = DB_postspherical_ICs(obj);

elseif strcmp(id, 'radial')
    p = [0;0;0]-1e-15;

elseif strcmp(id, 'eccentric')
    [r, p] = DB_eccentric_ICs(obj);
end

return
