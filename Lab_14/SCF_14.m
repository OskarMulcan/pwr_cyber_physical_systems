
%% Zadanie 2

function y = f(x)
    y = sqrt(x.^2 + 1);
end

syms x;  
figure; 

x_in = linspace(0, 100, 500); 

y_actual = f(x_in);
plot(x_in, y_actual, '-', 'LineWidth', 2); 
hold on; 

x_points = {2, 3}
for x_0 = x_points
    y_lin = taylor(f(x), x, 'ExpansionPoint', x_0, 'Order', 2); 
    y_lin_func = matlabFunction(y_lin);
    
    y_lin_out = y_lin_func(x_in); 
    
    plot(x_in, y_lin_out, '--');
end 

xlabel('x'); 
ylabel('y');
title('Multi-point Linearization of f(x) = \sqrt{x^2 + 1}');
grid on;
legend('Actual Function f(x)', 'Linear Approximations', 'Location', 'best');