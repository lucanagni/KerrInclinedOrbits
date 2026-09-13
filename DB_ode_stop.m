function [value,isterminal,direction] = DB_ode_stop(T,Y,rend)

%==========================================================================
% Stopping condition for the dynamics
% Assumes cartesian coordinates
% PR: 06/02/2019
%==========================================================================

r = sqrt(Y(1).^2 + Y(2).^2 + Y(3).^2);
value     = r-rend;
isterminal= 1;
direction = 0;

end

