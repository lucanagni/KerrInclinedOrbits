function [n,lamb,l] = DB_coprecessing_basis(dyn)
    stop = length(dyn.t);
    R = [dyn.x,dyn.y,dyn.z];

    n = zeros(stop,3);
    l = zeros(stop,3);
    lamb = zeros(stop,3);
    
    dndt = zeros(stop,3);

    for i=1:stop
        n(i,:) = R(i,:)./norm(R(i,:));
    end
    for i=1:3
        dndt(:,i) = DB_D1(n(:,i),dyn.t,4);
    end 
    for i=1:stop
        lamb(i,:) = dndt(i,:)./norm(dndt(i,:));
        l(i,:) = cross(n(i,:),lamb(i,:));
    end

return