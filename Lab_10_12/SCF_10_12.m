%{
Systemy Cyber-Fizyczne - Laboratorium 10-12
Zjawisko chaosu deterministycznego w modelowaniu systemów
%}

close all; clc; clear;

%% Plots
params.x0_vector = [0.01, 0.1, 0.3, 0.301, 0.5, 0.699, 0.7, 0.9, 0.99];
params.m = 30;
plot_logistic_equation(params.m, 0.9, params.x0_vector);
plot_logistic_equation(params.m, 1.5, params.x0_vector);
plot_logistic_equation(params.m, 2.8, params.x0_vector);
plot_logistic_equation(params.m, 3.4, params.x0_vector);
plot_logistic_equation(params.m, 3.6, params.x0_vector);

params.x0_step = 0.002;
params.x0_values = params.x0_step:params.x0_step:1-params.x0_step;
plot_le_bifurcation_diagram(0:0.005:4, params.x0_values, 200, 150);
plot_le_bifurcation_diagram(3.5:0.001:4, params.x0_values, 200, 150);

%% Functions
function [out] = logistic_equation(m, alpha, x0)
    arguments
        m     (1, 1) {mustBeNumeric, mustBeInteger, mustBeGreaterThanOrEqual(m, 1)}
        alpha (1, 1) {mustBeNumeric, mustBeBetween(alpha, 0, 4)}
        x0    (1, 1) {mustBeNumeric, mustBeBetween(x0, 0, 1)}
    end

    out = zeros(1, m + 1);
    out(1) = x0;
    
    for n = 1:m
        out(n + 1) = alpha * out(n) * (1 - out(n));
    end
end

function [] = plot_logistic_equation(m, alpha, x0_vector)
    arguments
        m         (1, 1) {mustBeNumeric, mustBeInteger, mustBeGreaterThanOrEqual(m, 1)}
        alpha     (1, 1) {mustBeNumeric, mustBeBetween(alpha, 0, 4)}
        x0_vector (1, :) {mustBeNumeric, mustBeBetween(x0_vector, 0, 1)}
    end

    figure;
    hold on;
    
    n           = length(x0_vector);
    generations = 0:m;
    plot_colors = hsv(n);

    for i = 1:n
        x0 = x0_vector(i);
        plot(generations, logistic_equation(m, alpha, x0), ...
             '-o', LineWidth=1.5, MarkerSize=2, ...
             DisplayName=['x_0 = ' num2str(x0)], ...
             Color=plot_colors(i, :));
    end

    xlabel('Generation (n)');
    ylabel('Population (x_n)');
    title(['Logistics equation (\alpha = ' num2str(alpha) ')']);

    ylim([0 1]);
    grid on;
    legend('show', 'Location', 'southoutside', 'Orientation', 'horizontal', 'NumColumns', 6);
    hold off;
end

function [] = plot_le_bifurcation_diagram(alpha_vector, x0_vector, m_stabilization, m_plot)
    arguments
        alpha_vector    (1, :) {mustBeNumeric, mustBeBetween(alpha_vector, 0, 4)}
        x0_vector       (1, :) {mustBeNumeric, mustBeBetween(x0_vector, 0, 1)}
        m_stabilization (1, 1) {mustBeNumeric, mustBeInteger, mustBeGreaterThanOrEqual(m_stabilization, 1)}
        m_plot          (1, 1) {mustBeNumeric, mustBeInteger, mustBeGreaterThanOrEqual(m_plot, 1)}
    end

    plot_color = [0, 0.28, 0.67];
    m_all = m_stabilization + m_plot;

    figure;
    hold on;

    for x0 = x0_vector
        all_alpha = zeros(1, length(alpha_vector) * m_plot);
        all_x     = zeros(1, length(alpha_vector) * m_plot);

        idx = 1;
        for alpha = alpha_vector
            out = logistic_equation(m_all, alpha, x0);
            last_gens = out((m_stabilization + 2):end);
            all_alpha(idx:(idx + m_plot - 1)) = alpha;
            all_x(idx:(idx + m_plot - 1))     = last_gens;
            idx = idx + m_plot;
        end

        scatter(all_alpha, all_x, 0.6, 'filled', ...
            'MarkerFaceColor', plot_color, ...
            'MarkerFaceAlpha', 0.075);
    end

    xlabel('Parameter \alpha');
    ylabel('Population (x_\infty)');
    title('Bifurcation Diagram for Logistic Equation');
    xlim([min(alpha_vector) max(alpha_vector)]);
    ylim([0 1]);
    grid on;

    hold off;
end

alpha_values = 0:0.005:4;
n_transient = 1000;
n_iterations = 500;
lambda = zeros(size(alpha_values));

for i = 1:length(alpha_values)
    a = alpha_values(i);
    x = 0.5;

    for j = 1:n_transient
        x = a * x * (1 - x);
    end
    
    sum_log = 0;
    for j = 1:n_iterations
        deriv = abs(a * (1 - 2 * x));
        sum_log = sum_log + log(deriv);
        
        x = a * x * (1 - x);
    end
    
    lambda(i) = sum_log / n_iterations;
end

figure;
plot(alpha_values, lambda, 'k', 'LineWidth', 0.5);
hold on;
yline(0, 'r--', 'LineWidth', 1.5);
xlabel('\alpha');
ylabel('\lambda');
title('Wykładnik Lapunowa dla odwzorowania logistycznego');
grid on;
ylim([-2 1]);