function f = freq(h,t)
    phi = -unwrap(angle(h));
    f = DB_D1(phi,t,4);
return