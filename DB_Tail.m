function t = DB_Tail(l,k,hatk,bphys,varargin)
%------------------------------------------------
% This function computes the tail contribution in
% the comparable mass case
%------------------------------------------------
% USAGE: EOBTail(l,k,hatk,bphys) where it is
%
%                l     -> the multipole
%                k     -> m Omega
%                hatk  -> m Omega Hreal
%                bphys -> 2
% @IHES 15/05/2007

ratio   = DB_GammaComplex(l+1-2.0*1i*hatk)./DB_GammaComplex(l+1);

if isempty(varargin)
    t = ratio.*exp(pi.*hatk).*exp(2.d0*1i*hatk.*log(2.d0*k.*bphys));
else
    alpha = varargin{1}
    t = ratio.*exp(pi.*hatk.*alpha).*exp(2.d0*1i*hatk.*log(2.d0*k.*bphys));
end

return
