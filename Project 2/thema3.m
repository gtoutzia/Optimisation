
syms x y;
f = x.^3 .* exp(-x.^2 - y.^4);

[x1,y1,fmin1,k1] = newton_gamma_const(f,-1 ,-1,0.01,0.1);
fprintf("(%4f,%4f) , minvalue=%4f after %d steps with constant gamma\n\n",x1(k1),y1(k1),fmin1,k1);
[x2,y2,fmin2,k2] = newton_gamma_min(f,-1,-1,0.01);
fprintf("(%4f,%4f) , minvalue=%4f after %d steps with min gamma\n\n",x2(k2),y2(k2),fmin2,k2);
[x3,y3,fmin3,k3] = newton_armijo(f,-1,-1,0.01,0.01,0.4);
fprintf("(%4f,%4f) , minvalue=%4f after %d steps with armijo gamma\n\n",x3(k3),y3(k3),fmin3,k3);


function [xk, yk, fmin, k] = newton_gamma_const(f, x0, y0, epsilon, gammak)
    syms x y;
    k = 1;
    xk = zeros(1, 1000);  
    yk = zeros(1, 1000);  
    xk(1) = x0;
    yk(1) = y0;
    f_values = zeros(1, 1000);
    fmin = NaN;
    gradf = gradient(f, [x, y]);
    hessf = hessian(f, [x, y]);
    
    while norm(double(subs(gradf, [x, y], [xk(k), yk(k)]))) >= epsilon
        gradient_at_point = double(subs(gradf, [x, y], [xk(k), yk(k)]));
        gradsqr_at_point = double(subs(hessf, [x, y], [xk(k), yk(k)]));
        inv_hess = inv(gradsqr_at_point);
        
        dkx = -inv_hess(1, 1) * gradient_at_point(1);
        dky = -inv_hess(2, 2) * gradient_at_point(2);
        
        
        xk(k + 1) = xk(k) + gammak * dkx;
        yk(k + 1) = yk(k) + gammak * dky;
        
        fmin = double(subs(f, [x, y], [xk(k), yk(k)]));
        f_values(k) = fmin;
        k = k+1;
    end
    
    figure;
    plot(1:k-1, f_values(1:k-1), 'o-');
    xlabel('Number of Steps (k)');
    ylabel('f Values');
    title('Convergence of f during Newton Optimization using constant step');
    grid on;
    end

function [xk, yk, fmin, k] = newton_gamma_min(f, x0, y0, epsilon)
    syms x y;
    k = 1;
    xk = zeros(1, 1000);  
    yk = zeros(1, 1000);  
    xk(1) = x0;
    yk(1) = y0;
    gradf = gradient(f, [x, y]);
    hessf = hessian(f, [x, y]);
    fmin = NaN;
    f_values = zeros(1, 1000);
    while norm(double(subs(gradf, [x, y], [xk(k), yk(k)]))) >= epsilon
        if(k > 50)
            break
        end
        
        gradient_at_point = double(subs(gradf, [x, y], [xk(k), yk(k)]));
        gradsqr_at_point = double(subs(hessf, [x, y], [xk(k), yk(k)]));
        
        inv_hess = inv(gradsqr_at_point);
        dk(1) = -inv_hess(1, 1) * gradient_at_point(1);
        dk(2) = -inv_hess(2, 2) * gradient_at_point(2);
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
    

     
        xk(k+1) = xk(k) + gammak * dk(1);  
        yk(k+1) = yk(k) + gammak * dk(2);

   
        fmin = double(subs(f, [x, y], [xk(k), yk(k)]));
        f_values(k) = fmin;
        
        k = k + 1;
    end

    figure;
    plot(1:k-1, f_values(1:k-1), 'o-');
    xlabel('Number of Steps (k)');
    ylabel('f Values');
    title('Convergence of f during Newton Optimization using best step');
    grid on;
    end



function [xk, yk, fmin, k] = newton_armijo(f, x0, y0, epsilon, alpha, beta)
    syms x y;
    k = 1;
    xk = zeros(1, 1000);
    yk = zeros(1, 1000);
    xk(1) = x0;
    yk(1) = y0;
    fmin = NaN;
    gradf = gradient(f, [x, y]);
    hessf = hessian(f, [x, y]);
    f_values = zeros(1, 1000);
    while norm(double(subs(gradf, [x, y], [xk(k), yk(k)]))) >= epsilon
        gradient_at_point = double(subs(gradf, [x, y], [xk(k), yk(k)]));
        gradsqr_at_point = double(subs(hessf, [x, y], [xk(k), yk(k)]));
        inv_hess = inv(gradsqr_at_point);

        dkx = -inv_hess(1, 1) * gradient_at_point(1);
        dky = -inv_hess(2, 2) * gradient_at_point(2);

        
        t = 1; 
        while double(subs(f, [x, y], [xk(k) + t * dkx, yk(k) + t * dky])) > ...
              double(subs(f, [x, y], [xk(k), yk(k)]) + alpha * t * dot(gradient_at_point, [dkx; dky]))
            t = beta * t;  
        end

        xk(k + 1) = xk(k) + t * dkx;
        yk(k + 1) = yk(k) + t * dky;
        fmin = double(subs(f, [x, y], [xk(k), yk(k)]));
        f_values(k) = fmin;
        k = k + 1;
        if k >100
            break
        end
    end  
    figure;
    plot(1:k-1, f_values(1:k-1), 'o-');
    xlabel('Number of Steps (k)');
    ylabel('f Values');
    title('Convergence of f during Newton Optimization using Armijo step');
    grid on;
    end

