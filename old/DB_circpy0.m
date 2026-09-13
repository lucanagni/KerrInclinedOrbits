function py0 = DB_circpy0(obj)

%==========================================================================
% Calculates initial quasispherical py0
% Takes r0,q,chi1,chi2 as input
% PR: 07/02/2019
%
% SA: updated with class logic 
%==========================================================================

r0    = obj.r0;
pphi0 = 6;
pphicirc = fzero(@(pphi) dHeff_dr(obj,r0,pphi),pphi0);

py0 = pphicirc./r0;

%{ 
%DEBUG -> plots dHeff_dr as a function of pphi (check zeros)

pyspan = linspace(0,10,100);
%py = linspace*0;
for i=1:length(pyspan)
    py(i) = dHeff_dr(obj,r0,pyspan(i));
end
figure
plot(pyspan,py);
hold on
yline(0);
title(obj.hamiltonian)
%}

return

function dHeff_dr = dHeff_dr(obj, x, pphi)

r = [x;0;0];
p = [0;pphi./x;0];

[~,~,dHeff] = DB_Hamiltonian(obj,r,p,obj.chi1,obj.chi2); %Hamiltonian is called with p instead of ps even if kerr_star because with these initial cond they are equivalent

if contains(obj.hamiltonian, 'star') %this shouldn't be needeed for the same reason
    [T,~] = DB_Tmatrix(r,obj.chi1);
    dHeff.dp = (dHeff.dps'*T)';
end

dHeff_dr = dHeff.dx(1,:) - 1./x.^2.*pphi.*dHeff.dp(2,:);

return