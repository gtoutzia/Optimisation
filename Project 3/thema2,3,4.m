clc;
syms x y;
f = 1/3 * x^2 + 3*y^2;

steepest_descent_projection(f,5,-5,-10,5,-8,12,5,0.01,0.5);
steepest_descent_projection(f,-5,10,-10,5,-8,12,15,0.01,0.1);
steepest_descent_projection(f,8,-10,-10,5,-8,12,0.1,0.01,0.2);
function [xk, yk, k] = steepest_descent_projection(f, x0, y0, a1, a2, b1, b2, sk, epsilon, gammak)
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

        if xk(k) - sk * gradient_at_point(1) <= a1
            x_bar = a1;
        elseif xk(k) - sk * gradient_at_point(1) >= a2
            x_bar = a2;
        else
            x_bar = xk(k) - sk * gradient_at_point(1);
        end

        if yk(k) - sk * gradient_at_point(2) <= b1
            y_bar = b1;
        elseif yk(k) - sk * gradient_at_point(2) >= b2
            y_bar = b2;
        else
            y_bar = yk(k) - sk * gradient_at_point(2);
        end

        % Update step with projection
        xk(k + 1) = max(a1, min(xk(k) - gammak * (double(subs(gradf(1), [x, y], [xk(k), yk(k)])) - x_bar), a2));
        yk(k + 1) = max(b1, min(yk(k) - gammak * (double(subs(gradf(2), [x, y], [xk(k), yk(k)])) - y_bar), b2));

        f_values(k) = double(subs(f, [x, y], [xk(k), yk(k)]));
        k = k + 1;
        if k > 100
            break
        end
        fmin = double(subs(f, [x, y], [xk(k), yk(k)]));
    end

    % Plot convergence
    figure;
    plot(1:k-1, f_values(1:k-1), 'o-');
    xlabel('Number of Steps (k)');
    ylabel('f Values');
    title(sprintf('Convergence of f during steepest descent Optimization using constant step gamma=%f', gammak));
    grid on;

    fprintf('Fmin=%f at (%f,%f)\n', fmin, xk(k), yk(k));
end



    