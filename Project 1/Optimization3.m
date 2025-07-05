
a1 = 0;
b1 = 3;
l_values1 = 0.01:0.01:0.1;
l_values2 = [0.01, 0.05, 0.1];

plot_steps_vs_l(@f1, a1, b1, l_values1);
plot_steps_vs_l(@f2, a1, b1, l_values1);
plot_steps_vs_l(@f3, a1, b1, l_values1);

plot_intervals_for_l(@f1, a1, b1, l_values2);
plot_intervals_for_l(@f2, a1, b1, l_values2);
plot_intervals_for_l(@f3, a1, b1, l_values2);

function [a,b,k] = fibonacci(f, a1, b1, l)
    fib = [1 1];
    while fib(end) < (b1-a1)/l
        fib(end+1) = fib(end) + fib(end-1);
    end
    k = 1;
    n = length(fib);
    x1 = a1 + (fib(n-2)/fib(n))*(b1-a1);
    x2 = a1 + (fib(n-1)/fib(n))*(b1-a1);
    a = [a1 zeros(1, n-2)];
    b = [b1 zeros(1, n-2)];
    for i = 1:n-2
        d = (fib(n-i-1)/fib(n-i))*(b(i)-a(i));
        if f(x1) < f(x2)
            b(i+1) = x2;
            x2 = x1;
            x1 = a(i+1) + d;
        else
            a(i+1) = x1;
            x1 = x2;
            x2 = b(i+1) - d;
        end
        k = i + 1;
    end
    a = a(1:k);
    b = b(1:k);
end

function plot_intervals_for_l(f, a1, b1, l_values)
    figure;
    for i = 1:length(l_values)
        l = l_values(i);
        [a, b, k] = fibonacci(f, a1, b1, l);
        
        subplot(length(l_values), 1, i);
        plot(1:k, a, '-o', 1:k, b, '-o');
        xlabel('k');
        ylabel('Values');
        functionName = func2str(f);
        title(['a and b for l = ' num2str(l) ' (' functionName ')']);
        legend('a', 'b');
        grid on;
    end
    set(gcf, 'Position', [100, 100, 800, 800]);
    sgtitle('Fibonacci Method Iterations');
end

function plot_steps_vs_l(f, a1, b1, l_values)
    steps = zeros(size(l_values));
    for i = 1:length(l_values)
        l = l_values(i);
        [~, ~, k] = fibonacci(f, a1, b1, l);
        steps(i) = k;
    end
    figure;
    plot(l_values, steps, '-o');
    xlabel('l');
    ylabel('Number of Steps');
    functionName = func2str(f);
    title(['Number of Steps vs. l (Fibonacci Method) - ' functionName]);
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
