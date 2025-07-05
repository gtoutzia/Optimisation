%thema1
clc;
syms x y;
f = x.^3 .* exp(-x.^2 - y.^4);
[x, y] = meshgrid(-2:0.1:2, -2:0.1:2);
f_numeric = matlabFunction(f);
z = f_numeric(x, y);
figure;
surf(x, y, z);
xlabel('X-axis');
ylabel('Y-axis');
zlabel('f(x, y)');
title('3D Surface Plot of f(x, y) = x^3 * exp(-x^2 - y^4)');
