function [t, x, v] = velocity_verlet(f, x_0, v_0, t_i, t_f, h)
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
        x(k+1) = x(k) + v(k)*h + (0.5)*(h^2)*f(t(k), x(k), v(k));

        v(k+1) = v(k) + h*0.5*(f(t(k), x(k), v(k)) + f(t(k+1), x(k+1), v(k+1)));
    end

end