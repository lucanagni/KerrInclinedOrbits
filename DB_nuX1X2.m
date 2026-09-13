function [nu, X1, X2] = DB_nuX1X2(q)
    X1 = q/(1+q);
    X2 = 1/(1+q);
    nu = X1*X2;
return

