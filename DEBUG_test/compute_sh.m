function [Y_teuk,Y_wig] = compute_sh(theta,phi,l)

    for i=1:(2*l+1)
        m = i-3;
        disp(l)
        disp(m)
        pause
        Y_teuk(i) = spinsphericalharm(0,l,m,phi,theta);
        Y_wig(i) = sphericalharm(l,m,phi,theta);
    end
return

% Definition of Ylm
function Y = spinsphericalharm(s, l, m, phi, i)
    if l<0 || m<-l || m>l
        error("wrong (l,m) inside spinsphericalharm(s,l,m,phi,i)");
    end
    c       = (-1).^(s)*sqrt( (2.*l+1.)./(4.*pi) );
    dWigner = c.*wigner_d_function(l,m,-s,i);
    rY      = cos(m*phi).*dWigner;
    iY      = sin(m*phi).*dWigner; %originally was + sign (see https://arxiv.org/pdf/0709.0093 eq (11.7) THAT IS THE CORRECT DEFINITION!
    Y       = rY + 1i*iY;
return
    
function out = wigner_d_function(l,m,s,i)  %as written, this computes d^l_{m,s} (order of indices is important)
    cth  = cos(i*0.5);
    sth  = sin(i*0.5);
    norm = sqrt( (factorial(l+m) * factorial(l-m) * factorial(l+s) * factorial(l-s)) );
    ki   = max(0,s-m);
    kf   = min(l+s,l-m);
    dWig = 0;
    for k=ki:kf
        div  = 1.0./( factorial(k) * factorial(l+s-k) * factorial(l-k-m) * factorial(k-s+m) );
        dWig = dWig + div*( (-1).^(k-s+m) * cth.^(2*l-m+s-2*k) * sth.^(2*k-s+m) );
    end
    out = dWig*norm;
return

function out = wigner_d_function2(l,m,s,i)  %as written, this computes d^l_{m,s} (order of indices is important)
    cth  = cos(i*0.5);
    sth  = sin(i*0.5);
    norm = sqrt( (factorial(l+m) * factorial(l-m) * factorial(l+s) * factorial(l-s)) );
    ki   = max(0,s-m);
    kf   = min(l+s,l-m);
    dWig = 0;
    for k=ki:kf
        div  = 1.0./( factorial(k) * factorial(l+s-k) * factorial(l-k-m) * factorial(k-s+m) );
        dWig = dWig + div*( (-1).^(k) * cth.^(2*l-m+s-2*k) * sth.^(2*k-s+m) );
    end
    out = dWig*norm;
return

function Y = sphericalharm(l, m, phi, i)
    if l<0 || m<-l || m>l
        error("wrong (l,m) inside spinsphericalharm(s,l,m,phi,i)");
    end
    c       = (-1).^(m)*sqrt( (2.*l+1.)./(4.*pi) );
    dWigner = c.*wigner_d_function2(l,m,0,i);
    rY      = cos(m*phi).*dWigner;
    iY      = sin(m*phi).*dWigner; %originally was + sign (see https://arxiv.org/pdf/0709.0093 eq (11.7) THAT IS THE CORRECT DEFINITION!
    Y       = rY + 1i*iY;
return