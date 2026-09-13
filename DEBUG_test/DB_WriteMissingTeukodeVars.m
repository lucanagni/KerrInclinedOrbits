function DB_WriteMissingTeukodeVars(t, rp, thp, php, ...
                                      drp, dthp, dphp, ...
                                      ddrp, ddthp, ddphp, ...
                                      dddrp, dddthp, dddphp, ...
                                      put_kin, dput_kin, ddput_kin, dddput_kin, ...
                                      pur_kin, dpur_kin, ddpur_kin, dddpur_kin, ...
                                      puth_kin, dputh_kin, ddputh_kin, dddputh_kin, ...
                                      puph_kin, dpuph_kin, ddpuph_kin, dddpuph_kin, ...
                                      pdt_kin, ...
                                      pdr_kin, ...
                                      pdth_kin, ...
                                      pdph_kin, ...
                                      mbh, abh, mu, spin, ...
                                      outputdirname, steps)                                    
%          
% REMINDER: previous function declaration, using down-index momenta
% --------------------------------------------------------------------
% function kerr_WriteMissingTeukodeVars(t, rp, thp, php, ...
%   drp, dthp, dphp, ...
%   ddrp, ddthp, ddphp, ...
%   dddrp, dddthp, dddphp, ...
%   Pt, dPt, ddPt, dddPt, ...
%   Pr, dPr, ddPr, dddPr,...
%   Pth, dPth, ddPth, dddPth, ...
%   Pph, dPph, ddPph, dddPph, ...
%   mbh, abh, mu, spin, ...
%   outputdirname, steps)
% --------------------------------------------------------------------
%
% NOTE: many things here may be redundant now, that we get UP-INDEX
% kinematical momenta passed ...
%
% this function computes from the vars that are produced by
% the KerrOrbitSolver, the vars that are needed in the Teukode, i.e.
% 
%  1-t
%  2-r     3-dr     4-ddr     5-dddr
%  6-th    7-dth    8-ddth    9-dddth
% 10-ph   11-dph   12-ddph   13-dddph
% 14-v^t  15-dv^t  16-ddv^t  17-dddv^t
% 18-v^r  19-dv^r  20-ddv^r  21-dddv^r
% 22-v^th 23-dv^th 24-ddv^th 25-dddv^th
% 26-v^ph 27-dv^ph 28-ddv^ph 29-dddv^ph
% 30-p^t  31-dp^t  32-ddp^t  33-dddp^t
% 34-p^r  35-dp^r  36-ddp^r  37-dddp^r
% 38-p^th 39-dp^th 40-ddp^th 41-dddp^th
% 42-p^ph 43-dp^ph 44-ddp^ph 45-dddp^ph
% 46-S^01 47-dS^01 48-ddS^01 49-dddS^01
% 50-S^02 51-dS^02 52-ddS^02 53-dddS^02
% 54-S^03 55-dS^03 56-ddS^03 57-dddS^03
% 58-S^12 59-dS^12 60-ddS^12 61-dddS^12
% 62-S^13 63-dS^13 64-ddS^13 65-dddS^13
% 66-S^23 67-dS^23 68-ddS^23 69-dddS^23
%                      


% further vars
abh2   = abh.*abh;
rp2    = rp.*rp;  
sinyp  = sin(thp);
cscyp  = 1./sinyp;
sinyp2 = sinyp.*sinyp;
cscyp2 = cscyp.*cscyp;
sin2yp = sin(2.*thp);
sin3yp = sin(3.*thp);
cosyp  = cos(thp); 
cos2yp = cos(2.*thp);
cos3yp = cos(3.*thp);
cotyp  = cosyp./sinyp; 
HRZp   = rp2 - 2.*mbh.*rp + abh2;
HRZp2  = HRZp.*HRZp;
HRZp3  = HRZp2.*HRZp;
HRZp4  = HRZp2.*HRZp2;
Sryp   = rp2 + abh2.*cosyp.*cosyp;
Sryp2  = Sryp.*Sryp;  
Sryp3  = Sryp2.*Sryp;
Sryp4  = Sryp2.*Sryp2;

% Renaming (to use my other formulas)
Vur = drp;
Vuy = dthp;
Vup = dphp;

dVur  = ddrp;
dVuy = ddthp;
dVup = ddphp;

ddVur  = dddrp;
ddVuy = dddthp;
ddVup = dddphp;

% powers
Vur2  = Vur.*Vur;
Vuy2  = Vuy.*Vuy;
drp2  = Vur2; % somehow few have not been replace to Vuyr2,Vuy2 correctly
dthp2 = Vuy2;

%=========================================================================%
%============================= Compute up index spin vector, tensor
%                                               and drvts 
%=========================================================================%

% for CEO we start from the assumption:
% 
%      S^mu = (0,0,S^th,0)
%
% ---> S_mu = (0,0,S_th,0)
%
%      Sigma^2 = -S_mu S^mu = -gdthdth S^th S^th 
%      
%  --> S^th =  S./sqrt(-gdthdth)
%      S_th = sqrt(gdthdth) .* S
% 
% then:   S^ab = epsilon^abcd S_c p_d./mass_p
%
% , where epsilon^abcd = epssym^abcd .* det(-g)^-1
%
% with epssym^0123 = -1 .
%
%  --> S^tr   = epsilon^0123 S_th p_ph = -det(-g)^-0.5 .* S_th .* p_ph
%      S^tth  = epsilon^02cd S_c  p_d  = 0
%      S^tph  = epsilon^0321 S_th p_r  = det(-g)^-0.5 .* S_th .* p_r
%      S^rth  = epsilon^12cd S_c  p_d  = 0
%      S^rph  = epsilon^1320 S_th p_t  = -det(-g)^-0.5 .* S_th .* p_t
%      S^thph = epsilon^23cd S_c  p_d  = 0

% spin-vector 
% Sut = zeros(steps,1);
% Sur = zeros(steps,1);
% Suy = mu*spin./Sryp;
% Sup = zeros(steps,1);

% drvts pf spin vector
% ACHTUNG: TODO will be switched on, for anything that is not a CEO
% dSut = zeros(steps,1);
% dSur = zeros(steps,1);
% dSuy = zeros(steps,1); 
% dSup = zeros(steps,1);
% 
% ddSut = zeros(steps,1);
% ddSur = zeros(steps,1);
% ddSuy = zeros(steps,1); 
% ddSup = zeros(steps,1);
% 
% dddSut = zeros(steps,1);
% dddSur = zeros(steps,1);
% dddSuy = zeros(steps,1); 
% dddSup = zeros(steps,1);

% aux down index spin vector
% gdthdth_p = -Sryp;
% Sdy = gdthdth_p.*Suy;
sqrtmdetgBL = sinyp.*Sryp;

% spin tensor (using -+++) metric
Suutr  = -spin.*pdph_kin.*sqrt(Sryp)./sqrtmdetgBL;
Suuty  = zeros(steps,1); 
Suutp  = spin.*pdr_kin.*sqrt(Sryp)./sqrtmdetgBL;
Suury  = zeros(steps,1); 
Suurp  = -spin.*pdt_kin.*sqrt(Sryp)./sqrtmdetgBL;
Suuyp  = zeros(steps,1); 

% drvts of spin tensor (all zero for CEO !)
% ACHTUNG: TODO will be switched on, for anything that is not a CEO
dSuutr  = zeros(steps,1); 
dSuuty  = zeros(steps,1); 
dSuutp  = zeros(steps,1); 
dSuury  = zeros(steps,1); 
dSuurp  = zeros(steps,1); 
dSuuyp  = zeros(steps,1); 

ddSuutr  = zeros(steps,1); 
ddSuuty  = zeros(steps,1); 
ddSuutp  = zeros(steps,1); 
ddSuury  = zeros(steps,1); 
ddSuurp  = zeros(steps,1); 
ddSuuyp  = zeros(steps,1); 

dddSuutr  = zeros(steps,1); 
dddSuuty  = zeros(steps,1); 
dddSuutp  = zeros(steps,1); 
dddSuury  = zeros(steps,1); 
dddSuurp  = zeros(steps,1); 
dddSuuyp  = zeros(steps,1); 

%
% ACHTUNG: now we get directly the up-index kinematical momenta. so it is
% not needed anymore to compute them here from the down-index ones.
% 
% ( do not delete, maybe needed later on again ..)
%
%=========================================================================%
%============================================== Compute up index momenta
%                                               and drvts 
%=========================================================================%
%
%
% 
% % upindex metric functions (+---)
% gutut_p    = (2.*abh2.*mbh.*rp.*sinyp2 + (abh2 + rp2).*Sryp)./(HRZp.*Sryp);
% gurur_p    = -HRZp./Sryp;
% guthuth_p  = -1./Sryp;
% guphuph_p  = -((cscyp2.*(HRZp - abh2.*sinyp2))./(HRZp.*Sryp)) ;
% gutuph_p   = (2.*abh.*mbh.*rp)./(HRZp.*Sryp);
% 
% % first compute up index momenta 
% put_kin        = gutut_p   .* Pt  + gutuph_p.*Pph;	
% pur_kin        = gurur_p   .* Pr  ;
% puth_kin       = guthuth_p .* Pth ;
% puph_kin       = guphuph_p .* Pph + gutuph_p.*Pt;
% 
% %%%%%%%%%%%%%%%%% 1st time drvts %%%%%%%%%%%%%%%%%%%%%%%%
% dHRZp      = 2.*Vur.*(-mbh + rp);
% dHRZp2     = dHRZp.*dHRZp;
% dHRZp3     = dHRZp2.*dHRZp;
% dSryp      = 2.*Vur.*rp - abh2.*Vuy.*sin2yp;
% dSryp2     = dSryp.*dSryp;
% dSryp3     = dSryp2.*dSryp;
% 
% dgutut_p     = (-(dHRZp.*Sryp.*(2.*abh2.*mbh.*rp.*sinyp2 + abh2.*Sryp + rp2.*Sryp)) + ...
%      2.*HRZp.*(-(abh2.*dSryp.*mbh.*rp.*sinyp2) + ...
%         abh2.*mbh.*sinyp.*(2.*cosyp.*Vuy.*rp + Vur.*sinyp).*Sryp + Vur.*rp.*Sryp2))./(HRZp2.*Sryp2);
% dgurur_p     = (dSryp.*HRZp - dHRZp.*Sryp)./Sryp2;
% dguthuth_p   = dSryp./Sryp2;
% dguphuph_p   = (dSryp.*HRZp.*(-abh2 + cscyp2.*HRZp) + (-(abh2.*dHRZp) + 2.*cotyp.*cscyp2.*Vuy.*HRZp2).*Sryp)./ ...
%    (HRZp2.*Sryp2);
% dgutuph_p    = (-2.*abh.*mbh.*(dHRZp.*rp.*Sryp + HRZp.*(dSryp.*rp - Vur.*Sryp)))./(HRZp2.*Sryp2);
% 
% dput_kin       = dPph.*gutuph_p  + dPt.*gutut_p + dgutuph_p.*Pph + dgutut_p.*Pt;   
% dpur_kin       = dgurur_p   .*Pr  + dPr.*gurur_p;   
% dputh_kin      = dguthuth_p .*Pth + dPth.*guthuth_p;
% dpuph_kin      = dPph.*guphuph_p + dPt.*gutuph_p + dguphuph_p.*Pph + dgutuph_p.*Pt;
% 
% %%%%%%%%%%%%%%%%% 2nd time drvts %%%%%%%%%%%%%%%%%%%%%%%%
% ddHRZp       = 2.*(Vur2 + dVur.*(-mbh + rp)); 
% ddSryp       = 2.*Vur2 - 2.*abh2.*cos2yp.*Vuy2 + 2.*dVur.*rp - abh2.*dVuy.*sin2yp; 
%   
% ddgutut_p    = (2.*dHRZp2.*(2.*abh2.*mbh.*rp.*sinyp2 + abh2.*Sryp + rp2.*Sryp).*Sryp2 - ...
%     HRZp.*Sryp.*(-4.*abh2.*dHRZp.*dSryp.*mbh.*rp.*sinyp2 + ...
%         2.*abh2.*mbh.*sinyp.*(2.*dHRZp.*Vur.*sinyp + ...
%            rp.*(4.*cosyp.*dHRZp.*Vuy + ddHRZp.*sinyp)).*Sryp + ...
%         (abh2.*ddHRZp + 4.*dHRZp.*Vur.*rp + ddHRZp.*rp2).*Sryp2) + ...
%      2.*HRZp2.*(2.*abh2.*dSryp2.*mbh.*rp.*sinyp2 - ...
%         abh2.*mbh.*sinyp.*(2.*Vur.*dSryp.*sinyp + ...
%            rp.*(4.*cosyp.*dSryp.*Vuy + ddSryp.*sinyp)).*Sryp + ...
%         abh2.*mbh.*(2.*Vur.*Vuy.*sin2yp + rp.*(2.*cos2yp.*Vuy2 + dVuy.*sin2yp) + ...
%            dVur.*sinyp2).*Sryp2 + (Vur2 + dVur.*rp).*Sryp3))./(HRZp3.*Sryp3);
% ddgurur_p    = (-2.*dSryp2.*HRZp + 2.*dHRZp.*dSryp.*Sryp + ddSryp.*HRZp.*Sryp - ddHRZp.*Sryp2)./Sryp3;
% ddguthuth_p  = (-2.*dSryp2 + ddSryp.*Sryp)./Sryp3;
% ddguphuph_p  = (2.*dSryp2.*(abh2 - cscyp2.*HRZp).*HRZp2 + ...
%      HRZp.*(2.*abh2.*dHRZp.*dSryp + ...
%         HRZp.*(-(abh2.*ddSryp) + cscyp2.*(ddSryp - 4.*cotyp.*dSryp.*Vuy).*HRZp)).* ...
%       Sryp + (2.*abh2.*dHRZp2 - HRZp.* ...
%          (abh2.*ddHRZp + (cscyp.^4).*HRZp2.* ...
%             (2.*(2 + cos2yp).*Vuy2 - dVuy.*sin2yp))).*Sryp2)./(HRZp3.*Sryp3);
% ddgutuph_p   = (2.*abh.*mbh.*(-(HRZp.*Sryp.*(-2.*dHRZp.*dSryp.*rp + (2.*dHRZp.*Vur + ddHRZp.*rp).*Sryp)) +  ...
%        2.*dHRZp2.*rp.*Sryp2 + HRZp2.*(2.*dSryp2.*rp - (2.*Vur.*dSryp + ddSryp.*rp).*Sryp + dVur.*Sryp2) ...
%        ))./(HRZp3.*Sryp3) ; 
%        
% ddput_kin      = 2.*dgutuph_p.*dPph + 2.*dgutut_p.*dPt + ddPph.*gutuph_p + ddPt.*gutut_p +  ...
%    ddgutuph_p.*Pph + ddgutut_p.*Pt;   
% ddpur_kin      = 2.*dgurur_p.*dPr + ddPr.*gurur_p + ddgurur_p.*Pr;   
% ddputh_kin     = 2.*dguthuth_p.*dPth + ddPth.*guthuth_p + ddguthuth_p.*Pth;
% ddpuph_kin     = 2.*dguphuph_p.*dPph + 2.*dgutuph_p.*dPt + ddPph.*guphuph_p +  ...
%    ddPt.*gutuph_p + ddguphuph_p.*Pph + ddgutuph_p.*Pt;
% 
% 
% %%%%%%%%%%%%%%%%% 3rd time drvts %%%%%%%%%%%%%%%%%%%%%%%%
% dddHRZp       = 2.*ddVur.*(-mbh + rp) + 6.*dVur.*Vur; 
% dddSryp       = 2.*ddVur.*rp - abh2.*ddVuy.*sin2yp + 6.*dVur.*Vur - 6.*abh2.*cos2yp.*dVuy.*Vuy + ...
%    4.*abh2.*sin2yp.*(Vuy.^3); 
%   
% dddgutut_p    = (-((2.*abh2.*mbh.*rp.*sinyp2 + abh2.*Sryp + rp2.*Sryp).* ...
%         (6.*dSryp3.*HRZp3 + 6.*dSryp.*(dHRZp.*dSryp - ddSryp.*HRZp).*HRZp2.*Sryp +  ...
%           HRZp.*(6.*dHRZp2.*dSryp - 3.*ddSryp.*dHRZp.*HRZp + ...
%              HRZp.*(-3.*ddHRZp.*dSryp + dddSryp.*HRZp)).*Sryp2 + ...
%           (6.*dHRZp3 - 6.*ddHRZp.*dHRZp.*HRZp + dddHRZp.*HRZp2).*Sryp3)) + ...
%      3.*HRZp.*Sryp.*(2.*dSryp2.*HRZp2 + HRZp.*(2.*dHRZp.*dSryp - ddSryp.*HRZp).*Sryp + ...
%         (2.*dHRZp2 - ddHRZp.*HRZp).*Sryp2).* ...
%       (abh2.*dSryp + dSryp.*rp2 + 2.*(abh2.*mbh.*sinyp2 + rp.*Sryp).*Vur + ...
%         2.*abh2.*mbh.*rp.*sin2yp.*Vuy) - 3.*HRZp2.*(dSryp.*HRZp + dHRZp.*Sryp).*Sryp2.* ...
%       (abh2.*ddSryp + abh2.*dVur.*mbh - abh2.*cos2yp.*dVur.*mbh + 4.*abh2.*cos2yp.*dthp2.*mbh.*rp + ...
%         ddSryp.*rp2 + 2.*abh2.*dVuy.*mbh.*rp.*sin2yp + 2.*(drp2 + dVur.*rp).*Sryp +  ...
%         4.*Vur.*(dSryp.*rp + abh2.*mbh.*sin2yp.*Vuy)) + ...
%      HRZp3.*Sryp3.*(6.*drp2.*dSryp + dddSryp.*rp2 + ...
%         6.*(2.*abh2.*cos2yp.*dthp2.*mbh + ddSryp.*rp + abh2.*dVuy.*mbh.*sin2yp + dVur.*Sryp).*Vur + ...
%         abh2.*(dddSryp + 2.*ddVur.*mbh.*sinyp2 + 6.*dVur.*mbh.*sin2yp.*Vuy) + ...
%         2.*rp.*(3.*dSryp.*dVur + abh2.*ddVuy.*mbh.*sin2yp + ddVur.*Sryp + ...
%            6.*abh2.*cos2yp.*dVuy.*mbh.*Vuy - 4.*abh2.*mbh.*sin2yp.*(Vuy.^3))))./(HRZp4.*Sryp4);
% dddgurur_p    = (6.*dSryp3.*HRZp - 6.*dSryp.*(dHRZp.*dSryp + ddSryp.*HRZp).*Sryp +  ...
%      (3.*ddSryp.*dHRZp + 3.*ddHRZp.*dSryp + dddSryp.*HRZp).*Sryp2 - dddHRZp.*Sryp3)./Sryp4;
% dddguthuth_p  = (6.*dSryp3 - 6.*ddSryp.*dSryp.*Sryp + dddSryp.*Sryp2)./Sryp4;
% dddguphuph_p  = (6.*dSryp3.*(-abh2 + cscyp2.*HRZp).*HRZp3 + ...
%      Sryp3.*(-6.*abh2.*dHRZp3 + 6.*abh2.*ddHRZp.*dHRZp.*HRZp - abh2.*dddHRZp.*HRZp2 + ...
%         ((cscyp.^5).*HRZp4.*(4.*cosyp.*ddVuy.*sinyp2 - 6.*dVuy.*(sin3yp + 3.*sinyp).*Vuy + ...
%              4.*(cos3yp + 11.*cosyp).*(Vuy.^3)))./2.) + ...
%      6.*dSryp.*HRZp2.*Sryp.*(-(abh2.*dHRZp.*dSryp) + ...
%         HRZp.*(abh2.*ddSryp + cscyp2.*HRZp.*(-ddSryp + 2.*cotyp.*dSryp.*Vuy))) + ...
%      HRZp.*Sryp2.*(-6.*abh2.*dHRZp2.*dSryp + 3.*abh2.*ddSryp.*dHRZp.*HRZp + ...
%         HRZp.*(3.*abh2.*ddHRZp.*dSryp - abh2.*dddSryp.*HRZp + ...
%            ((cscyp.^4).*HRZp2.*(12.*(2 + cos2yp).*dSryp.*dthp2 + ...
%                 2.*sinyp.*(-6.*cosyp.*dSryp.*dVuy + dddSryp.*sinyp) - 6.*ddSryp.*sin2yp.*Vuy))./2.)))./ ...
%    (HRZp4.*Sryp4);
% dddgutuph_p   = (2.*abh.*mbh.*(-6.*dHRZp3.*rp.*Sryp3 + 6.*dHRZp.*HRZp.*Sryp2.* ...
%         (-(dHRZp.*dSryp.*rp) + Sryp.*(ddHRZp.*rp + dHRZp.*Vur)) +  ...
%        HRZp3.*(-6.*dSryp3.*rp + ddVur.*Sryp3 - ...
%           Sryp2.*(3.*dSryp.*dVur + dddSryp.*rp + 3.*ddSryp.*Vur) + ...
%           6.*dSryp.*Sryp.*(ddSryp.*rp + dSryp.*Vur)) - ...
%        HRZp2.*Sryp.*(6.*dHRZp.*dSryp2.*rp + Sryp2.*(3.*dHRZp.*dVur + dddHRZp.*rp + 3.*ddHRZp.*Vur) - ...
%           3.*Sryp.*(ddSryp.*dHRZp.*rp + ddHRZp.*dSryp.*rp + 2.*dHRZp.*dSryp.*Vur))))./(HRZp4.*Sryp4); 
%  
% dddput_kin      = 3.*ddPph.*dgutuph_p + 3.*ddPt.*dgutut_p + 3.*ddgutuph_p.*dPph + ...
%    3.*ddgutut_p.*dPt + dddPph.*gutuph_p + dddPt.*gutut_p + ...
%    dddgutuph_p.*Pph + dddgutut_p.*Pt;   
% dddpur_kin      = 3.*ddPr.*dgurur_p + 3.*ddgurur_p.*dPr + dddPr.*gurur_p + ...
%    dddgurur_p.*Pr;   
% dddputh_kin     = 3.*ddPth.*dguthuth_p + 3.*ddguthuth_p.*dPth + ...
%    dddPth.*guthuth_p + dddguthuth_p.*Pth;
% dddpuph_kin     = 3.*ddPph.*dguphuph_p + 3.*ddPt.*dgutuph_p + ...
%    3.*ddguphuph_p.*dPph + 3.*ddgutuph_p.*dPt + ...
%    dddPph.*guphuph_p + dddPt.*gutuph_p + dddguphuph_p.*Pph + ...
%    dddgutuph_p.*Pt;
   
   

%=========================================================================%
%============================================== Compute up index 4-velocity
%                                               and drvts 
%=========================================================================%


%%% NOTE: when building p^mu from p_mu the signature of the metric is
%%% important!
%%% --> not done anymore, so commmented this sign correction here
%
%%% have used (+---) metric above as in orig Teuk paper and initially in
%%% teukode
%%% --> seems better to use (-+++) 
% put_kin      = -put_kin;   
% pur_kin      = -pur_kin;   
% puth_kin     = -puth_kin;
% puph_kin     = -puph_kin;
% 
% dput_kin     = -dput_kin;   
% dpur_kin     = -dpur_kin;   
% dputh_kin    = -dputh_kin;
% dpuph_kin    = -dpuph_kin;
% 
% ddput_kin    = -ddput_kin;   
% ddpur_kin    = -ddpur_kin;   
% ddputh_kin   = -ddputh_kin;
% ddpuph_kin   = -ddpuph_kin;
% 
% dddput_kin   = -dddput_kin;   
% dddpur_kin   = -dddpur_kin;   
% dddputh_kin  = -dddputh_kin;
% dddpuph_kin  = -dddpuph_kin;


%%% Linear in spin approximation (cf. Eq.(2.10) in [Faye, Blanchet, Buonanno,2006] )
vut      = put_kin;   
vur      = pur_kin;   
vuth     = puth_kin;
vuph     = puph_kin;

dvut     = dput_kin;   
dvur     = dpur_kin;   
dvuth    = dputh_kin;
dvuph    = dpuph_kin;

ddvut    = ddput_kin;   
ddvur    = ddpur_kin;   
ddvuth   = ddputh_kin;
ddvuph   = ddpuph_kin;

dddvut   = dddput_kin;   
dddvur   = dddpur_kin;   
dddvuth  = dddputh_kin;
dddvuph  = dddpuph_kin;


%=========================================================================%
%==============================================        Output  
%=========================================================================%

% readable format

%fid(1)  = fopen([outputdirname,'/vut.dat']      ,'w');
%fid(2)  = fopen([outputdirname,'/vur.dat']     ,'w');
%fid(3)  = fopen([outputdirname,'/vur.dat']     ,'w');
%fid(4)  = fopen([outputdirname,'/vuph.dat']    ,'w');
fid(5)  = fopen([outputdirname,'/put_kin.dat']     ,'w');
fid(6)  = fopen([outputdirname,'/pur_kin.dat']     ,'w');
fid(7)  = fopen([outputdirname,'/pur_kin.dat']     ,'w');
fid(8)  = fopen([outputdirname,'/puph_kin.dat']    ,'w');
fid(9)  = fopen([outputdirname,'/Str.dat']     ,'w');
fid(10) = fopen([outputdirname,'/Stth.dat']    ,'w');
fid(11) = fopen([outputdirname,'/Stph.dat']    ,'w');
fid(12) = fopen([outputdirname,'/Srth.dat']    ,'w');
fid(13) = fopen([outputdirname,'/Srph.dat']    ,'w');
fid(14) = fopen([outputdirname,'/Sthph.dat']   ,'w');

for n=1:steps
    %fprintf(fid(1), '%15.12f %15.12f %15.12f  %15.12f %15.12f \n',[t(n) vut(n)  dvut(n)  ddvut(n)  dddvut(n)  ] ); 
    %fprintf(fid(2), '%15.12f %15.12f %15.12f  %15.12f %15.12f \n',[t(n) vur(n)  dvur(n)  ddvur(n)  dddvur(n)  ] );
    %fprintf(fid(3), '%15.12f %15.12f %15.12f  %15.12f %15.12f \n',[t(n) vuth(n) dvuth(n) ddvuth(n) dddvuth(n) ] );
    %fprintf(fid(4), '%15.12f %15.12f %15.12f  %15.12f %15.12f \n',[t(n) vuph(n) dvuph(n) ddvuph(n) dddvuph(n) ] );
    fprintf( fid(5), '%15.12f %15.12f %15.12f  %15.12f %15.12f \n',[t(n) put_kin(n)    dput_kin(n)    ddput_kin(n)    dddput_kin(n)  ] ); 
    fprintf( fid(6), '%15.12f %15.12f %15.12f  %15.12f %15.12f \n',[t(n) pur_kin(n)    dpur_kin(n)    ddpur_kin(n)    dddpur_kin(n)  ] );
    fprintf( fid(7), '%15.12f %15.12f %15.12f  %15.12f %15.12f \n',[t(n) puth_kin(n)   dputh_kin(n)   ddputh_kin(n)   dddputh_kin(n) ] );
    fprintf( fid(8), '%15.12f %15.12f %15.12f  %15.12f %15.12f \n',[t(n) puph_kin(n)   dpuph_kin(n)   ddpuph_kin(n)   dddpuph_kin(n) ] );
    fprintf( fid(9), '%15.12f %15.12f %15.12f  %15.12f %15.12f \n',[t(n) Suutr(n)  dSuutr(n)  ddSuutr(n)  dddSuutr(n)  ] ); 
    fprintf(fid(10), '%15.12f %15.12f %15.12f  %15.12f %15.12f \n',[t(n) Suuty(n)  dSuuty(n)  ddSuuty(n)  dddSuuty(n)  ] ); 
    fprintf(fid(11), '%15.12f %15.12f %15.12f  %15.12f %15.12f \n',[t(n) Suutp(n)  dSuutp(n)  ddSuutp(n)  dddSuutp(n)  ] ); 
    fprintf(fid(12), '%15.12f %15.12f %15.12f  %15.12f %15.12f \n',[t(n) Suury(n)  dSuury(n)  ddSuury(n)  dddSuury(n)  ] ); 
    fprintf(fid(13), '%15.12f %15.12f %15.12f  %15.12f %15.12f \n',[t(n) Suurp(n)  dSuurp(n)  ddSuurp(n)  dddSuurp(n)  ] ); 
    fprintf(fid(14), '%15.12f %15.12f %15.12f  %15.12f %15.12f \n',[t(n) Suuyp(n)  dSuuyp(n)  ddSuuyp(n)  dddSuuyp(n)  ] ); 
end

for k=5:14
    fclose(fid(k));
end

% One Column as read in by Teukode
% 
%  1-t
%  2-r     3-dr     4-ddr     5-dddr
%  6-th    7-dth    8-ddth    9-dddth
% 10-ph   11-dph   12-ddph   13-dddph
% 14-v^t  15-dv^t  16-ddv^t  17-dddv^t
% 18-v^r  19-dv^r  20-ddv^r  21-dddv^r
% 22-v^th 23-dv^th 24-ddv^th 25-dddv^th
% 26-v^ph 27-dv^ph 28-ddv^ph 29-dddv^ph
% 30-p^t  31-dp^t  32-ddp^t  33-dddp^t
% 34-p^r  35-dp^r  36-ddp^r  37-dddp^r
% 38-p^th 39-dp^th 40-ddp^th 41-dddp^th
% 42-p^ph 43-dp^ph 44-ddp^ph 45-dddp^ph
% 46-S^01 47-dS^01 48-ddS^01 49-dddS^01
% 50-S^02 51-dS^02 52-ddS^02 53-dddS^02
% 54-S^03 55-dS^03 56-ddS^03 57-dddS^03
% 58-S^12 59-dS^12 60-ddS^12 61-dddS^12
% 62-S^13 63-dS^13 64-ddS^13 65-dddS^13
% 66-S^23 67-dS^23 68-ddS^23 69-dddS^23
%  
nvars = 69;

OneColumnData = [ t; 
                  rp ; drp ; ddrp; dddrp; ...
                  thp ; dthp ; ddthp; dddthp; ...
                  php ; dphp ; ddphp; dddphp; ...
                  vut ; dvut ; ddvut; dddvut; ...
                  vur ; dvur ; ddvur; dddvur; ...
                  vuth; dvuth ; ddvuth; dddvuth; ...
                  vuph; dvuph ; ddvuph; dddvuph; ...
                  put_kin ; dput_kin ; ddput_kin; dddput_kin; ...
                  pur_kin ; dpur_kin ; ddpur_kin; dddpur_kin; ...
                  puth_kin; dputh_kin ; ddputh_kin; dddputh_kin; ...
                  puph_kin; dpuph_kin ; ddpuph_kin; dddpuph_kin; ...
                  Suutr; dSuutr; ddSuutr; dddSuutr; ...
                  Suuty; dSuuty; ddSuuty; dddSuuty; ...
                  Suutp; dSuutp; ddSuutp; dddSuutp; ...
                  Suury; dSuury; ddSuury; dddSuury; ...
                  Suurp; dSuurp; ddSuurp; dddSuurp; ...
                  Suuyp; dSuuyp; ddSuuyp; dddSuuyp; ];                  
if(length(OneColumnData) ~= nvars*steps )
  disp(' Dimensions not right.. sth. wrong in building data for teukode.');
  return;
end
 
% Write data
a_s  = sprintf('%.4f',abh);
r_s  = sprintf('%.2f',rp(1));
mu_s = sprintf('%.3e',mu);
Si_s = sprintf('%.2f',spin);

disp('Write long 1d data for teukode');

%fid(1)  = fopen([outputdirname,'/EOB_DN0mu',mu_s,'_a',a_s,'_S',Si_s,'_r',r_s,'.dat'],'w');
fid(1)  = fopen([outputdirname,'/SP_CEO_BDHAM_a',a_s,'_S',Si_s,'_r',r_s,'_mu',mu_s,'.dat'],'w');
fprintf(fid(1),'#nt %i\n',steps);
fprintf(fid(1),'#nv 69\n');
fprintf(fid(1),'#abh %.6f\n',abh);
fprintf(fid(1),'%.15e\n',OneColumnData);
fclose(fid(1));

fid(1)  = fopen(['./testshitdata.dat'],'w');
fprintf(fid(1),'%.15e\n',[1,2,3;1,2,3;3,4,5]);
fclose(fid(1));


end
