%{
Systemy Cyber-Fizyczne - Laboratorium 6 i 7
Modelowanie ruchu falowego
%}

close all; clc; clear;

%% Zadanie 1 - obliczenia

function [h] = h_xt_forloop(x, t, c)
    arguments
        x (1, :) {mustBeNumeric}
        t (1, :) {mustBeNumeric}
        c (1, 1) {mustBeNumeric}
    end
    L  = x(end);
    dt = t(2) - t(1); % zakładamy const dt
    dx = x(2) - x(1); % zakładamy const dx
    r  = c * dt / dx;
    nx = length(x);
    nt = length(t);
    h  = zeros(nx, nt);

    % Warunek początkowy t=0
    for x_idx = 1:nx
        h(x_idx, 1) = 0.5 - 0.5 * cos(2 * pi / L * x(x_idx));
    end

    % Warunek początkowy t=dt
    %{
        Zakładamy zerową prędkość początkową => dh/dt(0) = 0
        Przybliżamy pochodną różnicą centralną: dh/dt ~=~ (h(x, t+dt) - h(x, t-dt))/(2*dt)
        => dh/dt(0) = (h(x, dt) - h(x, -dt))/(2*dt) = 0
        => h(x, dt) = h(x, -dt)
        => h(x, dt) = r^2 * (h(x+dx, 0) + h(x-dx, 0)) + 2 * (1 - r^2) * h(x, 0) - h(x, dt)
        => h(x, dt) = (r^2 / 2) * (h(x+dx, 0) + h(x-dx, 0)) + (1 - r^2) * h(x, 0)
    %}
    for x_idx = 2:nx - 1
        h(x_idx, 2) = (r^2 / 2) * (h(x_idx + 1, 1) + h(x_idx - 1, 1)) + (1 - r^2) * h(x_idx, 1);
    end

    for t_idx = 2:nt - 1
        for x_idx = 2:nx - 1
            h(x_idx, t_idx + 1) = r^2 * (h(x_idx + 1, t_idx) + h(x_idx - 1, t_idx)) + 2 * (1 - r^2) * h(x_idx, t_idx) - h(x_idx, t_idx - 1);
        end
    end
end

zad1.x = 0:0.1:10;
zad1.t = 0:0.005:30;
zad1.c = 5;

zad1.h_forloop = h_xt_forloop(zad1.x, zad1.t, zad1.c);

zad1.points = zeros(length(zad1.t), length(zad1.x), 2);
zad1.points(:, :, 1) = repmat(zad1.x, length(zad1.t), 1); % x stałe dla każdego t
zad1.points(:, :, 2) = zad1.h_forloop';

%% Zadanie 1 - animacja
animate_2d(zad1.t, zad1.points, 'x', 'h(x,t)', 'Struna gitarowa');

%% Funkcje
function animate_2d(t, points, x_label, y_label, title_str)
    arguments
        t         (1, :)    {mustBeNumeric}
        points    (:, :, 2) {mustBeNumeric}  % (t_idx, point_idx, [x, y])
        x_label   (1, :)    char = 'x'
        y_label   (1, :)    char = 'y'
        title_str (1, :)    char = ''
    end

    dt = t(2) - t(1);

    figure;
    hLine = plot(squeeze(points(1, :, 1)), squeeze(points(1, :, 2)), 'b-', 'LineWidth', 2);
    xlim([min(points(:, :, 1), [], 'all'), max(points(:, :, 1), [], 'all')]);
    ylim([min(points(:, :, 2), [], 'all'), max(points(:, :, 2), [], 'all')]);
    xlabel(x_label); ylabel(y_label);
    grid on;

    for t_idx = 1:length(t)
        tic;
        hLine.XData = squeeze(points(t_idx, :, 1));
        hLine.YData = squeeze(points(t_idx, :, 2));
        if isempty(title_str)
            title(sprintf('t = %.2f s', t(t_idx)));
        else
            title(sprintf('%s, t = %.2f s', title_str, t(t_idx)));
        end
        drawnow;

        elapsed = toc;
        remaining = dt - elapsed;
        if remaining > 0
            pause(remaining);
        end
    end
end
