clc;
clf;
syms x y;
f = 1/3 * x^2 + 3*y^2;
[x, y] = meshgrid(-2:0.1:2, -2:0.1:2);
f_numeric = matlabFunction(f);
z = f_numeric(x, y);
figure;
surf(x, y, z);
xlabel('X-axis');
ylabel('Y-axis');
zlabel('f(x, y)');
title('3D Surface Plot of f');

steepest_descent(f,1,1,0.001,0.1);
steepest_descent(f,1,1,0.001,0.3);
steepest_descent(f,1,1,0.001,3);
steepest_descent(f,1,1,0.001,5);


function [xk, yk, fmin, k] = steepest_descent(f, x0, y0, epsilon, gammak)
    syms x y;
    k = 1;
    xk = zeros(1, 1000);  
    yk = zeros(1, 1000);  
    xk(1) = x0;
    yk(1) = y0;
    f_values = zeros(1, 1000);
    gradf = gradient(f, [x, y]);
    
    while norm(double(subs(gradf, [x, y], [xk(k), yk(k)]))) >= epsilon
        gradient_at_point = double(subs(gradf, [x, y], [xk(k), yk(k)]));
        xk(k + 1) = xk(k) - gammak * gradient_at_point(1);
        yk(k + 1) = yk(k) - gammak * gradient_at_point(2);
        fmin= double(subs(f, [x, y], [xk(k), yk(k)]));
        f_values(k) = fmin;
        k = k + 1;
        if k > 50 
            break
        end
    end
      fmin= double(subs(f, [x, y], [xk(k), yk(k)]));
    figure;
    plot(1:k-1, f_values(1:k-1), 'o-');
    xlabel('Number of Steps (k)');
    ylabel('f Values');
    title(sprintf('Convergence of f during steepest decent Optimization using constant step gamma=%f',gammak));
    grid on;
    fprintf('Fmin=%f at (%f,%f)\n',fmin,xk(k),yk(k));
end