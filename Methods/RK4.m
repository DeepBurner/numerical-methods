function [u, t] = RK4(f, u0, time, h)
    t0 = time(1);
    tn = time(2);
    
    n = ceil((tn-t0)/h)+1;
    
    t = zeros(n, 1);
    u = zeros(n, length(u0));
    
    t(1) = t0;
    u(1,:) = u0;
    
    for k = 1:(n-1)
        t(k+1) = t(k) + h;
        
        k1 = f(t(k), u(k,:));
        k2 = f(t(k) + 1/2*h, u(k,:) + 1/2*k1*h);
        k3 = f(t(k) + 1/2*h, u(k,:) + 1/2*k2*h);
        k4 = f(t(k) + h, u(k,:) + k3*h);
        u(k+1,:) = u(k,:) + h/6*(k1 + 2*k2 + 2*k3 + k4);
    end;
 end