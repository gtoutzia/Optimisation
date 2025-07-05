
l = 0.01;
epsilon = 0.001;
a1 = 0;
b1 = 3;
l_values= 0.001:0.01:0.01;

[af1,bf1,kf1]=bisection_method(@(x) f1(x),epsilon,l,a1,b1);
[af2,bf2,kf2]=bisection_method(@(x) f2(x),epsilon,l,a1,b1);
[af3,bf3,kf3]=bisection_method(@(x) f3(x),epsilon,l,a1,b1);

k_values1 = 1:kf1-1; 
k_values2 =1:kf2-1;
k_values3 =1:kf3-1;
figure;
subplot(3,1,1);
plot(k_values1, af1(1:kf1-1), '-o', k_values1, bf1(1:kf1-1), '-o');
xlabel('k');
ylabel('Values');
title('Bisection Method for f1(x)');
legend('af1', 'bf1');

subplot(3,1,2);
plot(k_values2, af2(1:kf2-1), '-o', k_values2, bf2(1:kf2-1), '-o');
xlabel('k');
ylabel('Values');
title('Bisection Method for f2(x)');
legend('af2', 'bf2');

subplot(3,1,3);
plot(k_values3, af3(1:kf3-1), '-o', k_values3, bf3(1:kf3-1), '-o');
xlabel('k');
ylabel('Values');
title('Bisection Method for f3(x)');
legend('af3', 'bf3');

plot_iterations(@(x) f1(x), 'f1(x)');
plot_iterations(@(x) f2(x), 'f2(x)');
plot_iterations(@(x) f3(x), 'f3(x)');
 
function [a,b,k] =bisection_method(f,epsilon,l,a1,b1)
    a(1) = a1;
    b(1) = b1;
    k = 1;
    while b(k) - a(k) >= l
        x1(k) = (a(k) + b(k)) / 2 - epsilon;
        x2(k) = (a(k) + b(k)) / 2 + epsilon;
       
        if f(x1(k)) < f(x2(k))
            a(k+1) = a(k);
            b(k+1) = x2(k);
        else 
            a(k+1) = x1(k);
            b(k+1) = b(k);
        end
        k = k + 1;
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

function k = count_iterations(f, epsilon, l, a1, b1)
    a = a1;
    b = b1;
    k = 0;
    while b - a >= l && k<15
        x1 = (a + b) / 2 - epsilon;
        x2 = (a + b) / 2 + epsilon;
        
        f_x1 = f(x1);
        f_x2 = f(x2);
        
        if f_x1 < f_x2
            a = a;
            b = x2;
        else 
            a = x1;
            b = b;
        end
        k = k + 1;
    end
end

function plot_iterations(f, function_name)
    epsilon_values = 0.001:0.001:0.01;
    l_values = 0.01:0.01:0.1;

    iterations_epsilon = zeros(size(epsilon_values));
    iterations_l = zeros(size(l_values));

    for i = 1:length(epsilon_values)
        iterations_epsilon(i) = count_iterations(f, epsilon_values(i), 0.01, 0, 3);
    end

    for i = 1:length(l_values)
        iterations_l(i) = count_iterations(f, 0.001, l_values(i), 0, 3);
    end

    % Plot number of iterations vs. epsilon
    figure;
    subplot(2,1,1);
    plot(epsilon_values, iterations_epsilon, '-o');
    xlabel('epsilon');
    ylabel('Number of Iterations');
    title(['Number of Iterations vs. epsilon for ' function_name]);
    grid on;

    % Plot number of iterations vs. l
    subplot(2,1,2);
    plot(l_values, iterations_l, '-o');
    xlabel('lambda');
    ylabel('Number of Iterations');
    title(['Number of Iterations vs. lambda for ' function_name]);
    grid on;
end