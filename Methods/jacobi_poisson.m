function [U, X, Y, h] = jacobi_poisson(init, final, N, bc, X_q, Y_q, Q, max_iter)
    
    %initialize constants & arrays
    e_0 = 8.854*10^(-12);

    x = linspace(init,final,N);
    y = linspace(init,final,N);
    h = x(2)-x(1);
    
    n = length(x);
    iter = 0;
    
    b = zeros(n,n);
    [X, Y] = meshgrid(x,y);
    U = zeros(n,n);
    
    %create the b matrix, loop over the space and put in b values wherever
    %applicable
    
    for i = 1:n-1
        for j = 1:n-1
            for k = 1:length(Q)
                if x(i) < X_q(k) && X_q(k) < x(i+1) && y(j) < Y_q(k) && Y_q(k) < y(j+1)

                    ii = ceil(2*(X_q(k)-x(i))/(x(i+1)-x(i))-1) + i;
                    jj = ceil(2*(Y_q(k)-y(j))/(y(i+1)-y(i))-1) + j;
                    
                    b(ii,jj) = -Q(k)/(e_0*h^2);
                end
            end;
        end;
    end;
    
    %boundary conditions
    U(:,1) = bc(1);
    U(1,:) = bc(2);
    U(:,end) = bc(3);
    U(end,:) = bc(4);

    %corner values
    U(1,1) = (bc(1)+bc(2))/2;
    U(1,end) = (bc(2)+bc(3))/2;
    U(end,1) = (bc(2)+bc(4))/2;
    U(end,end) = (bc(3)+bc(4))/2;
    
    %iterating over the solutions
    while iter < max_iter
        phi_updated = zeros(size(U));
        
        %update the values with the weighted averages
        for i = 2:(n-1)
            for j = 2:(n-1)
                phi_updated(i,j) = 0.25*( U(i+1,j) + U(i-1,j) +U(i,j+1) + U(i,j-1) - h^2*b(i,j));
            end;
        end;

        U = phi_updated;
        iter = iter + 1;
    end;

    U = U';
end