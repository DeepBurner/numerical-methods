function [t, x, v] = leapfrog(f, x_0, v_0, t_i, t_f, h)
    %%% Input arguments
    % f = second derivative of x with as a function of t, x and v
    % x_0 = initial value of x
    % v_0 = initial vlaue of v
    % t_i = starting time
    % t_f = finishing time
    % h = time step

    t = t_i:h:t_f;
    x = zeros(length(t));
    v = zeros(length(t));
    x(1) = x_0;
    v(1) = v_0;
    
    for k = 1:length(t)
        if k > 1
            v(k+1) = v(k-1) + 2*h*f(t(k), x(k), v(k));
            x(k+2) = x(k) + 2*h*v(k+1);
        else
            %not self starting
            v(k+1) = v(k) + h*f(t(k), x(k), v(k));
            x(k+1) = x(k) + h*v(k+1);
        end
    end

end