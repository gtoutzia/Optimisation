
a1 = 0;
b1 = 3;
l_values1 = 0.01:0.01:0.1;
l_values2 = [0.01, 0.05, 0.1];


df1 = @(x) (3*(x-1)^2 - 2*(x-4)*sin(x) - (x-4)^2*cos(x));
df2 = @(x) (-2*exp(-2*x) - 2*(x-2));
df3 = @(x) (2*x*log(0.5*x) + x - sin(0.2*x)^2 - 0.2*cos(0.2*x)^2);

plot_derivative_bisection_steps(@f1, df1, a1, b1, l_values1);
plot_derivative_bisection_steps(@f2, df2, a1, b1, l_values1);
plot_derivative_bisection_steps(@f3, df3, a1, b1, l_values1);

derivative_bisection_iterations_intervals(@f1, df1, a1, b1, l_values2);
derivative_bisection_iterations_intervals(@f2, df2, a1, b1, l_values2);
derivative_bisection_iterations_intervals(@f3, df3, a1, b1, l_values2);

function [a, b, k] = derivative_bisection(f, df, a1, b1, l)
    max_iterations = 1000;
    a = zeros(1, max_iterations);
    b = zeros(1, max_iterations);
    a(1) = a1;
    b(1) = b1;
    n = log(l / (b1 - a1)) / log(1/2);
    for k = 1:n
        x = (a(k) + b(k)) / 2;
        if df(x) == 0
            break;
        elseif df(x) > 0
            a(k+1) = a(k);
            b(k+1) = x;
        else
            a(k+1) = x;
            b(k+1) = b(k);
        end
    end
    a = a(1:k); 
    b = b(1:k);
end


function derivative_bisection_iterations_intervals(f, df, a1, b1, l_values)
    figure;
    for i = 1:length(l_values)
        l = l_values(i);
        [a, b, k] = derivative_bisection(f, df, a1, b1, l);
        
        subplot(length(l_values), 1, i);
        plot(1:k, a, 1:k, b);
        xlabel('k');
        ylabel('Values');
        functionName = func2str(f);
        title(['a and b for l = ' num2str(l) ' (' functionName ')']);
        legend('a', 'b');
        grid on;
    end
    set(gcf, 'Position', [100, 100, 800, 800]);
    sgtitle('Derivative Bisection Method Iterations');
end

function plot_derivative_bisection_steps(f, df, a1, b1, l_values)
    steps = zeros(size(l_values));
    
    for i = 1:length(l_values)
        l = l_values(i);
        [~, ~, k] = derivative_bisection(f, df, a1, b1, l);
        steps(i) = k;
    end
    
    figure;
    plot(l_values, steps, '-o');
    xlabel('l');
    ylabel('Number of Steps');
    functionName = func2str(f); 
    title(['Number of Steps vs. l (Derivative Bisection Method) - ' functionName]);
    grid on;
end

function y = f1(x)
    y = (x-1)^3 + (x-4)^2*cos(x);
end

function y = f2(x)
    y = exp(-2*x) + (x-2)^2;
end

function y = f3(x)
    y = x^2 * log(0.5*x) + sin(0.2*x)^2;
end


