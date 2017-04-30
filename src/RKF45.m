function [x,y,t] = RKF45(f,g,x0,y0,h0,steps, tol)
    %initialize the empty arrays

    t = zeros(1,steps);
    x = zeros(1,steps);
    y = zeros(1,steps);
    
    x(1) = x0;
    y(1) = y0;
    h = h0;
    
    %minimum and maximum values for timestep
    hmin = 0.000001;
    hmax = 0.01;
    
    k = 1;
    current_t = 0;
    
    while (k <= steps) %&& (current_t < 10)
        
        % check for hmax hmin
        if h < hmin
            h = hmin;
        elseif h > hmax
            h = hmax;
        end
        
        %RK45 coefficients
        
        kn1 = h * f(t(k), x(k), y(k));
        ln1 = h * g(t(k), x(k), y(k));
        
        kn2 = h * f(t(k) + h/4, x(k) + kn1/4, y(k) + ln1/4);
        ln2 = h * g(t(k) + h/4, x(k) + kn1/4, y(k) + ln1/4);
        
        kn3 = h * f(t(k) + 3*h/8, x(k) + 3*kn1/32 + 9*kn2/32, y(k) + 3*ln1/32 + 9*ln2/32);
        ln3 = h * g(t(k) + 3*h/8, x(k) + 3*kn1/32 + 9*kn2/32, y(k) + 3*ln1/32 + 9*ln2/32);
        
        kn4 = h * f(t(k) + 12*h/13, x(k) + 1932*kn1/2197 - 7200*kn2/2197 + 7296*kn3/2197, y(k) + 1932*ln1/2197 - 7200*ln2/2197 + 7296*ln3/2197);
        ln4 = h * g(t(k) + 12*h/13, x(k) + 1932*kn1/2197 - 7200*kn2/2197 + 7296*kn3/2197, y(k) + 1932*ln1/2197 - 7200*ln2/2197 + 7296*ln3/2197);
        
        kn5 = h * f(t(k) + h, x(k) + 439*kn1/216 - 8*kn2 + 3680*kn3/513 - 845*kn4/4104, y(k) + 439*ln1/216 - 8*ln2 + 3680*ln3/513 - 845*ln4/4104);
        ln5 = h * g(t(k) + h, x(k) + 439*kn1/216 - 8*kn2 + 3680*kn3/513 - 845*kn4/4104, y(k) + 439*ln1/216 - 8*ln2 + 3680*ln3/513 - 845*ln4/4104);
        
        kn6 = h * f(t(k) + h/2, x(k) - 8*kn1/27 + 2*kn2 - 3544*kn3/2565 + 1859*kn4/4104 - 11*kn5/40, y(k) - 8*ln1/27 + 2*ln2 - 3544*ln3/2565 + 1859*ln4/4104 - 11*ln5/40);
        ln6 = h * g(t(k) + h/2, x(k) - 8*kn1/27 + 2*kn2 - 3544*kn3/2565 + 1859*kn4/4104 - 11*kn5/40, y(k) - 8*ln1/27 + 2*ln2 - 3544*ln3/2565 + 1859*ln4/4104 - 11*ln5/40);
        
        xk_1 = x(k) + 25*kn1/216 + 1408*kn3/2565 + 2197*kn4/4101 - kn5/5;
        xzk_1 = x(k) + 16*kn1/135 + 6656*kn3/12825 + 28561*kn4/56430 - 9*kn5/50 + 2*kn6/55;
        
        yk_1 = y(k) + 25*ln1/216 + 1408*ln3/2565 + 2197*ln4/4101 - ln5/5;
        yzk_1 = y(k) + 16*ln1/135 + 6656*ln3/12825 + 28561*ln4/56430 - 9*ln5/50 + 2*ln6/55;
        
        %error amount
        error = abs(xzk_1 - xk_1);
        
        %check for error, adjust the timestep
        if error > tol
            h = h/2;
        else
            %x & y always >= 0
            if xk_1 < 0
                xk_1 = 0;
            end
            
            if yk_1 < 0
                yk_1 = 0;
            end
            
            t(k+1) = t(k) + h;
            x(k+1) = xk_1;
            y(k+1) = yk_1;
            
            
            current_t = current_t + h;
            k = k + 1; %only advance the sequence if conditions are met
            if error < tol/100
                h = 2 * h; %if it's too fine increase the time step
            end
        end
    end;
 end