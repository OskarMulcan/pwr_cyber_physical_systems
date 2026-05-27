%{
Systemy Cyber-Fizyczne - Laboratorium 10-12
Zjawisko chaosu deterministycznego w modelowaniu systemów
%}

close all; clc; clear;

%% Wykresy populacji
params.x0_vector = [0.01, 0.05, 0.1:0.2:0.9, 0.95, 0.99];
plot_logistic_equation(20, 1, params.x0_vector);
plot_logistic_equation(20, 1.5, params.x0_vector);
plot_logistic_equation(20, 2.8, params.x0_vector);
plot_logistic_equation(20, 3.4, params.x0_vector);
plot_logistic_equation(20, 3.6, params.x0_vector);

%% Funkcje
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

    generations = 0:m;

    for x0 = x0_vector
        plot(generations, logistic_equation(m, alpha, x0), ...
             '-o', LineWidth=1.5, MarkerSize=2, DisplayName=['x_0 = ' num2str(x0)]);
    end

    xlabel('Pokolenie (n)');
    ylabel('Populacja (x_n)');
    title(['Równanie logistyczne (\alpha = ' num2str(alpha) ')']);

    ylim([0 1]);
    grid on;
    legend('show', 'Location', 'best');
    hold off;
end