%thema1
clc;
syms x y;
f = x.^3 .* exp(-x.^2 - y.^4);

[x1,y1,fmin1,k1] = levenberg_constant_step(f, -1 , -1,0.1,0.1);
fprintf("Minimum was found using constant step at (x,y)= (%4f,%4f)\n with value %f after %d steps\n\n",x1(k1),y1(k1),fmin1,k1);
[x2,y2,fmin2,k2] = levenberg_bisection_step(f,-1,-1,0.1);
fprintf("Minimum was found using min step at (x,y)= (%4f,%4f)\n with value %f after %d steps\n\n",x2(k2),y2(k2),fmin2,k2);
[x3,y3,fmin3,k3] = levenberg_armijo_step(f,-1,-1,0.01,0.01,0.4);
fprintf("Minimum was found using armijo step at (x,y)= (%4f,%4f)\n with value %f after %d steps\n\n",x3(k3),y3(k3),fmin3,k3);

function [xk, yk, fmin, k] = levenberg_constant_step(f, x0, y0, epsilon, gammak)
    syms x y;
    k = 1;
    xk(1) = x0;
    yk(1) = y0;
    while true
        
        gradf = double(subs(gradient(f, [x, y]), [x, y], [xk(k), yk(k)]));
        hessf = double(subs(hessian(f , [x, y]), [x, y], [xk(k), yk(k)]));
        fval = double(subs(f, [x, y], [xk(k), yk(k)]));
        

        m = fval + 0.5;  
       I = eye(length(gradf));
        dk = linsolve(hessf + m*I, -gradf);

       
        xk(k + 1) = xk(k) + gammak * dk(1);
        yk(k + 1) = yk(k) + gammak * dk(2);
        fmin= double(subs(f, [x, y], [xk(k), yk(k)]));
        f_values(k) = fmin;
        if norm(gradf) < epsilon
            break;
        end

        k = k + 1;
    end

    fmin = double(subs(f, [x, y], [xk(k), yk(k)]));
   
    figure;
    plot(1:k-1, f_values(1:k-1), 'o-');
    xlabel('Number of Steps (k)');
    ylabel('f Values');
    title('Convergence of f during Levbemberg Optimization using constant step');
    grid on;
end

function [xk, yk, fmin, k] = levenberg_bisection_step(f, x0, y0, epsilon)
    syms x y;
    k = 1;
    xk(1) = x0;
    yk(1) = y0;

    while true
        gradf = double(subs(gradient(f, [x, y]), [x, y], [xk(k), yk(k)]));
        hessf = double(subs(hessian(f, [x, y]), [x, y], [xk(k), yk(k)]));
        fval = double(subs(f, [x, y], [xk(k), yk(k)]));

        m = fval + 0.5;
        I = eye(length(gradf));
        dk = linsolve(hessf + m*I, -gradf);

        
        phi = @(gamma) double(subs(f, [x, y], [xk(k) + gamma*dk(1), yk(k) + gamma*dk(2)]));

       
        a = 0;
        b = 1;
        tolerance = 1e-6;
        while (b - a) > tolerance
            gamma = (a + b) / 2;
            if phi(gamma - tolerance) < phi(gamma + tolerance)
                b = gamma;
            else
                a = gamma;
            end
        end

        gammak = (a + b) / 2;

        
        xk(k + 1) = xk(k) + gammak * dk(1);
        yk(k + 1) = yk(k) + gammak * dk(2);
        fmin= double(subs(f, [x, y], [xk(k), yk(k)]));
        f_values(k) = fmin;
        if norm(gradf) < epsilon
            break;
        end

        k = k + 1;
    end
    fmin = double(subs(f, [x, y], [xk(k), yk(k)]));
     
    figure;
    plot(1:k-1, f_values(1:k-1), 'o-');
    xlabel('Number of Steps (k)');
    ylabel('f Values');
    title('Convergence of f during Levemberg Optimization using bestg step');
    grid on;
end

function [xk, yk, fmin, k] = levenberg_armijo_step(f, x0, y0, epsilon, alpha, beta)
    syms x y;
    k = 1;
    xk(1) = x0;
    yk(1) = y0;
    gammak = 0.5;
    while true
        gradf = double(subs(gradient(f, [x, y]), [x, y], [xk(k), yk(k)]));
        hessf = double(subs(hessian(f, [x, y]), [x, y], [xk(k), yk(k)]));
        fval = double(subs(f, [x, y], [xk(k), yk(k)]));

        m = fval + 0.5;
        I = eye(length(gradf));
        dk = linsolve(hessf + m*I, -gradf);
       

       
        while double(subs(f, [x, y], [xk(k) + gammak*dk(1), yk(k) + gammak*dk(2)])) > fval + alpha*gammak*gradf'*dk
            gammak = beta * gammak;
        end

       
        xk(k + 1) = xk(k) + gammak * dk(1);
        yk(k + 1) = yk(k) + gammak * dk(2);
        fmin= double(subs(f, [x, y], [xk(k), yk(k)]));
        f_values(k) = fmin;
       
        if norm(gradf) < epsilon
            break;
        end

        k = k + 1;
    end
    fmin = double(subs(f, [x, y], [xk(k), yk(k)]));
    figure;
    plot(1:k-1, f_values(1:k-1), 'o-');
    xlabel('Number of Steps (k)');
    ylabel('f Values');
    title('Convergence of f during Levemberg Optimization using Armijo step');
    grid on;
end

