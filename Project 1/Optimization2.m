
a1 = 0;
b1 = 3;
l = 0.01;

l_values1 = 0.01:0.01:0.1;

plot_golden_ratio_steps(@f1, a1, b1, l_values1);
plot_golden_ratio_steps(@f2, a1, b1, l_values1);
plot_golden_ratio_steps(@f3, a1, b1, l_values1);

l_values2 = [0.01, 0.05, 0.1]; % Specify the tolerance levels
golden_ratio_iterations_intervals(@f1, a1, b1, l_values2);
golden_ratio_iterations_intervals(@f2, a1, b1, l_values2);
golden_ratio_iterations_intervals(@f3, a1, b1, l_values2);




function [a,b,k] = golden_ratio(f,a1,b1,l)
    gamma = 0.618;
    k =1;
    x1(k) = a1 + (1-gamma)*(b1 -a1);
    x2(k) = a1 +gamma*(b1 -a1);
    a(k) = a1;
    b(k) = b1;
    while b(k) - a(k) >= l
        if f(x1(k)) < f(x2(k))
            a(k+1) = x1(k);
            b(k+1) = b(k);
            x2(k+1) = a(k+1) +gamma*(b(k+1) - a(k+1));
            x1(k+1) =x2(k);
            k=k+1;
        
        else 
            a(k+1) = a(k);
            b(k+1) = x2(k);
            x2(k+1) = x1(k);
            x1(k+1) = a(k+1) +(1-gamma)*(b(k+1) -a(k+1));
            k = k+1;
        end
    
    end
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


function golden_ratio_iterations_intervals(f, a1, b1, l_values)
    gamma = 0.618;
    figure;
    
    for i = 1:length(l_values)
        l = l_values(i);
        [a, b, k] = golden_ratio(f, a1, b1, l);
        
        subplot(length(l_values), 1, i);
        plot(1:k, a(1:k), '-o', 1:k, b(1:k), '-o');
        xlabel('k');
        ylabel('Values');
        functionName = func2str(f);
        title(['a and b for l = ' num2str(l) ' (' functionName ')']);
        legend('a', 'b');
        grid on;
    end
    
    % Adjust figure properties
    set(gcf, 'Position', [100, 100, 800, 800]);
    sgtitle('Golden Ratio Method Iterations');
end

function plot_golden_ratio_steps(f, a1, b1, l_values)
    steps = zeros(size(l_values));
    
    for i = 1:length(l_values)
        l = l_values(i);
        [~, ~, k] = golden_ratio(f, a1, b1, l);
        steps(i) = k; % Store the step count for this l value
    end
    
    figure;
    plot(l_values, steps, '-o');
    xlabel('l');
    ylabel('Number of Steps');
    functionName = func2str(f); % Get the function name as a string
    title(['Number of Steps vs. l (Golden Ratio Method) - ' functionName]);
    grid on;
end



