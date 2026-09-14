function a = DB_testmass_checks(obj,varargin)

%==========================================================================
% Perform some tests when considering Kerr
% Return Kerr spin
%
% SA: 2/31/2024
%==========================================================================

if ~isempty(varargin)
    msg = char(varargin{1});
else
    msg = '';
end

q  = obj.q;
nu = q/(1+q)^2;

if nu>0.01
    error('Symmetric mass ratio too high for kerr dynamics: nu=%.5f\n', nu)
end

spin_kerr = obj.chi1;

if abs(spin_kerr(1))>1e-15 || abs(spin_kerr(2))>1e-15
    error('%sx and y component of Kerr spin are not zero!', msg)
end

a = spin_kerr(3);

return