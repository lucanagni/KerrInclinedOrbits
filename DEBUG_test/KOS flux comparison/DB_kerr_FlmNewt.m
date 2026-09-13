function s = DB_kerr_FlmNewt(x)
% DINFlmNewt.m 
%            This function computes the GW energy flux from a circularized
%            binary in its multipolar decoposition.
%            The output is given by Flm/nu^2, where nu=m1m2/(m1+m2)^2 is
%            the symmetric mass ratio. 
%            USAGE:
%            
%            F=DINFlmNewt(x,nu)
%
%            where F is a structure of the form s.lm. Values up to l=8 are
%            explicitly included.
%
%            Author: Alessandro Nagar
%            When:   @IHES, November 2008 & November 2009
%

% shorthands: variable
x5  = x.^5;
x6  = x.*x5;
x7  = x.*x6;
x8  = x.*x7;
x9  = x.*x8;
x10 = x.*x9;
x11 = x.*x10;
x12 = x.*x11;


% Newtonian partial fluxes

  
    s.l2m2    = x5* 32.d0/5.d0;                       %l2m2
    s.l2m1    = x6* 8.d0/45.d0;                       %l2m1
    s.l3m3    = x6* 243.d0/28.d0;                     %l3m3
    s.l3m2    = x7* 32.d0/63.d0;                      %l3m2
    s.l3m1    = x6/ 1260.d0;                          %l3m1
    s.l4m4    = x7* 8192.d0/567.d0;                   %l4m4
    s.l4m3    = x8* 729.d0/700.d0;                    %l4m3
    s.l4m2    = x7* 32.d0/3969.d0;                    %l4m2
    s.l4m1    = x8/ 44100.d0;                         %l4m1
    s.l5m5    = x8* 1953125.d0/76032.d0;              %l5m5
    s.l5m4    = x9* 131072.d0/66825.d0;               %l5m4
    s.l5m3    = x8* 2187.d0/70400.d0;                 %l5m3
    s.l5m2    = x9* 256.d0/400950.d0;                 %l5m2
    s.l5m1    = x8/ 19958400.d0;                      %l5m1
    s.l6m6    = x9* 839808.d0/17875.d0;               %l6m6
    s.l6m5    = x10* 48828125.d0/13621608.d0;         %l6m5
    s.l6m4    = x9* 4194304.d0/47779875.d0;           %l6m4
    s.l6m3    = x10* 59049.d0/15415400.d0;            %l6m3
    s.l6m2    = x9* 128.d0/28667925.d0;               %l6m2
    s.l6m1    = x10/ 1123782660.d0;                   %l6m1
    s.l7m7    = x10* 96889010407.d0/1111968000.d0;    %l7m7
    s.l7m6    = x11* 5668704.d0/875875.d0;            %l7m6
    s.l7m5    = x10* 1220703125.d0/5.666588928d9;     %l7m5
    s.l7m4    = x11* 4194304.d0/3.07432125d8;         %l7m4
    s.l7m3    = x10* 1594323.d0/3.2064032d10;         %l7m3
    s.l7m2    = x11* 32.d0/1.35270135d8;              %l7m2
    s.l7m1    = x10/ 9.3498717312d11;                 %l7m1
    s.l8m8    = x11* 274877906944.d0/1688511825.d0;   %l8m8
    s.l8m7    = x12* 4747561509943.d0/4.083146496d11; %l8m7
    s.l8m6    = x11* 51018336.d0/1.04229125d8;        %l8m6
    s.l8m5    = x12* 30517578125.d0/8.00296713216d11; %l8m5
    s.l8m4    = x11* 4194304.d0/1.5679038375d10;      %l8m4
    s.l8m3    = x12* 177147.d0/3.96428032d10;         %l8m3
    s.l8m2    = x11* 32.d0/3.4493884425d10;           %l8m2
    s.l8m1    = x12/ 8.174459284992d13;               %l8m0

return