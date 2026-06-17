%{
Systemy Cyber-Fizyczne - Laboratorium 14
Problematyka przybliżeń numerycznych i linearyzacji
%}

close all; clc; clear;

%% Zadanie 1
G = tf(5,[6 1 1]);

t.start = 0;
t.stop = 100;

syms t_sym s;
[num, den] = tfdata(G, 'v');
G_sym = poly2sym(num, s) / poly2sym(den, s);
y_t = ilaplace(G_sym / s, s, t_sym);
analytical_int = double(int(y_t, t_sym, t.start, t.stop));
fprintf('Analitycal: %.6f\n', analytical_int);

t.step_sizes = [1 0.5 0.1 0.05 0.001];

for step_size = t.step_sizes
    figure;
    [y, tout] = step(G, 0:step_size:100);

    euler_fwd_int = sum(y) * step_size;
    euler_bck_int = sum(y(1:end-1)) * step_size;
    trapez_int = sum((y(1:end-1) + y(2:end))) * step_size / 2;

    plot(tout, y, 'LineWidth', 1.5, 'DisplayName', 'Response y(t)');
    xlim([t.start, t.stop]);

    title(['Step Response for Step Size: ', num2str(step_size)]);
    xlabel('Time (s)');
    ylabel('Response');
    grid on;

    box_text = {
        ['Euler fwd: ', num2str(euler_fwd_int, '%.6f')], ...
        ['Euler bck: ', num2str(euler_bck_int, '%.6f')], ...
        ['Trapezoidal: ', num2str(trapez_int, '%.6f')]
        };

    annotation('textbox', [0.65, 0.15, 0.22, 0.15], ...
        'String', box_text, ...
        'BackgroundColor', [1 1 1 0.8], ...
        'EdgeColor', [0.5 0.5 0.5], ...
        'LineWidth', 1, ...
        'FontName', 'Helvetica', ...
        'FontSize', 10, ...
        'FitBoxToText', 'on');

    legend('Location','southeast');
end

%% Zadanie 2

function y = f(x)
    y = sqrt(x.^2 + 1);
end

syms x;  
figure; 

x_in = linspace(-100, 100, 500); 

y_actual = f(x_in);
plot(x_in, y_actual, '-', 'LineWidth', 2); 
hold on; 

x_points = {2, 3};
for x_0 = x_points
    y_lin = taylor(f(x), x, 'ExpansionPoint', x_0, 'Order', 2); 
    y_lin_func = matlabFunction(y_lin);
    
    y_lin_out = y_lin_func(x_in); 22
    
    plot(x_in, y_lin_out, '--');
end 

xlabel('x'); 
ylabel('y');
title('Multi-point Linearization of f(x) = \sqrt{x^2 + 1}');
grid on;
legend('Actual Function f(x)', 'Linear Approximations', 'Location', 'best');