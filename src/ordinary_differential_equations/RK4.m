function [t, x] = RK4(f, x_0, t_i, t_f, h)
    %%% Input arguments
    % f = second derivative of x with as a function of t, x and v
    % x_0 = initial value of x
    % t_i = starting time
    % t_f = finishing time
    % h = time step

    t = t_i:h:t_f;
    x = zeros(length(t));
    x(1) = x_0;
    
    for k = 1:length(t)
        
        k1 = h*f(t(k), x(k));
        k2 = h*f(t(k) + 1/2*h, x(k) + k1/2);
        k3 = h*f(t(k) + 1/2*h, x(k) + k2/2);
        k4 = h*f(t(k) + h, x(k) + k3);
        
        x(k+1) = x(k) + (k1 + 2*k2 + 2*k3 + k4)/6;
    end

end