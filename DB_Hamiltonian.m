function [Heff,Horb,dHeff,dHorb,dHso] = DB_Hamiltonian(obj, x, p, chi1, chi2)

%==========================================================================
% Wrapper for Hamiltonians
% SA: 01/31/24
%==========================================================================

q = obj.q;

if strcmp(obj.hamiltonian, 'balmelli')
    [Heff,H,dHeff,dH] = DB_Hamiltonian_Balmelli(x,p,q,chi1,chi2);

elseif strcmp(obj.hamiltonian, 'kerr')
    [Heff,Horb,dHeff,dHorb,dHso] = DB_Hamiltonian_Kerr(obj,x,p,chi1);

elseif strcmp(obj.hamiltonian, 'kerr_star')
    [Heff,Horb,dHeff,dHorb] = DB_Hamiltonian_Kerr_star(x,p,chi1);

elseif strcmp(obj.hamiltonian, 'kerr_eq')
    a = DB_testmass_checks(obj);
    [Heff,H,dHeff,dH] = DB_Hamiltonian_Kerr_eq(x,p,a);

else
    error('Unknown Hamiltonian: %s\n', obj.hamiltonian)
end

return
