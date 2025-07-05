
clc;
syms x y;
f = x.^3 .* exp(-x.^2 - y.^4);
[x1,y1,fmin1,k1] =max_steep_gamma_const(f,1,1,0.1,0.1);
fprintf("Minimum was found using constant step at (x,y)= (%4f,%4f)\n with value %f after %d steps\n\n",x1(k1),y1(k1),fmin1,k1);
[x2,y2,fmin2,k2] =max_steep_gamma_min(f,1,1,0.01);
fprintf("Minimum  was found using min step at (x,y)= (%4f,%4f)\n with value %f after %d steps\n\n",x2(k2),y2(k2),fmin2,k2);
[x3,y3,fmin3,k3] =max_steep_armijo(f,1,1,0.01,0.01,0.4);
fprintf("Minimum  was found using armijo step at (x,y)= (%4f,%4f)\n with value %f after %d steps\n\n",x3(k3),y3(k3),fmin3,k3);

function [xk, yk, fmin, k] = max_steep_gamma_const(f, x0, y0, epsilon, gammak)
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
    title('Convergence of f during steepest decent Optimization using constant step');
    grid on;
end


function [xk, yk,fmin,k] = max_steep_gamma_min(f,x0,y0,epsilon)
    syms x y;
    k = 1;
    xk = zeros(1, 1000);  
    yk = zeros(1, 1000); 
    xk(1) = x0;
    yk(1) = y0;
    f_values = zeros(1, 1000);
    gradf = gradient(f, [x, y]);
    gammak = 0.5;
    while norm(double(subs(gradf, [x, y], [xk(k), yk(k)]))) >= epsilon
        gradient_at_point = double(subs(gradf, [x, y], [xk(k), yk(k)]));
        gammak = bisection(xk(k) - gammak * gradient_at_point, 0.01,0.01,0,1);
        
        xk(k + 1) = xk(k) - gammak * gradient_at_point(1);
        yk(k + 1) = yk(k) - gammak * dgradient_at_point(2);
        fmin = double(subs(f, [x, y], [xk(k), yk(k)]));
        f_values(k)=fmin;
        k = k + 1;
        if k>50
            break
        end
    end
    fmin = double(subs(f, [x, y], [xk(k), yk(k)]));
    figure;
    plot(1:k-1, f_values(1:k-1), 'o-');
    xlabel('Number of Steps (k)');
    ylabel('f Values');
    title('Convergence of f during steepest decent Optimization using best step');
    grid on;
    end
    


function [xk, yk, fmin, k] = max_steep_armijo(f, x0, y0, epsilon, alpha, beta)
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
        step_size = 1.0;
        while subs(f, [x, y], [xk(k) - step_size * gradient_at_point(1), yk(k) - step_size * gradient_at_point(2)]) > ...
              subs(f, [x, y], [xk(k), yk(k)]) - alpha * step_size * dot(gradient_at_point, gradient_at_point)
            step_size = beta * step_size;
        end

        xk(k + 1) = xk(k) - step_size * gradient_at_point(1);
        yk(k + 1) = yk(k) - step_size * gradient_at_point(2);
        fmin = double(subs(f, [x, y], [xk(k), yk(k)]));
        f_values(k)=fmin;
        k = k + 1;
        if k>50
            break
        end
    end
    fmin = double(subs(f, [x, y], [xk(k), yk(k)]));
    figure;
    plot(1:k-1, f_values(1:k-1), 'o-');
    xlabel('Number of Steps (k)');
    ylabel('f Values');
    title('Convergence of f during steepest decent Optimization using Armijo step');
    grid on;
    
end


function min = bisection(f , e , l , a1 , b1)
k = 1;
ak(1,1) = a1;
bk(1,1) = b1;
if ((bk-ak)<l)
    return
else
    while((bk(k)-ak(k))>=l)
        x1k = (ak(k) + bk(k))/2 - e;
        x2k = (ak(k) + bk(k))/2 + e;
        if (subs(f,x1k)<subs(f,x2k))
            bk(k+1) = x2k;
            ak(k+1) = ak(k);
        else
            bk(k+1) = bk(k);
            ak(k+1) = x1k;
        end
        k = k + 1;
    end    
end    
min = (ak+bk)/2;
end
