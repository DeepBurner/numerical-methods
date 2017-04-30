%Ahmet Burak Catli
%11.3.2017
%PHYS371 Assignment 2: Polynomial Curve Fitting

%initialize the x and the J(x)
x = linspace(0, 3, 300);
J_x = inline('1 - 2.2499997.*(x./3).^2 + 1.2656208.*(x./3).^4 - 0.3163866.*(x/3).^6 + 0.0444479.*(x/3).^8 - 0.0039444.*(x./3).^10 + 0.0002100.*(x./3).^12', 'x');

%this table will be used to calculate the maximum deviation
table = zeros(11, 2);


%from polynomials of order 2 to 12
for m=2:12;
   
   %create the X matrix
   X = ones(300,1);
   for n=1:m
      X = [X, transpose(x).^n]; 
   end
   
   %create the y matrix
   y = transpose(J_x(x));
   C = transpose(X)*y;
   A = transpose(X)*X;
   
   %calculate theta(coefficient) matrix
   theta = linsolve(A,C);
   
   %calculate the different powers of x and add them up, fun is our final
   %approximation
   fun = zeros(1, 300);
   for z=0:m
      fun = fun + theta(z+1).*(x.^z);
      
      %plot the different powers
      plot(x, theta(z+1).*(x.^z));
      hold on;
   end
   
   %Find the maximum amount of deviation from original function
   table(m-1,1) = m;
   table(m-1,2) = max(abs(J_x(x) - fun));
   
   %Plot the final approximation function
   label = ['Approximation of J_x with max order of ' int2str(m)];
   plot(x, fun, '-.');
   hold on;
   title(label);
   xlabel('J_x');
   ylabel('x');
   
   %This command is for outputting the figures
   print(int2str(m), '-depsc2');
   clf;
end 