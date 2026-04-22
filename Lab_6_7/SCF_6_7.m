%{
Systemy Cyber-Fizyczne - Laboratorium 6 i 7
Modelowanie ruchu falowego
%}
close all; clc; clear;

%% Zadanie 1 - Struna gitarowa

% Parametry zadania 1
zad1.x = 0:0.1:10;
zad1.t = 0:0.005:30;
zad1.c = 5;
zad1.L = zad1.x(end);

% Warunki brzegowe h(0, t) i h(L, t):
zad1.bc_left = @(time) 0;
zad1.bc_right = @(time) 0;

% Warunek początkowy h(x, 0):
zad1.initial = @(x_val) 0.5 - 0.5 * cos(2 * pi / zad1.L * x_val);

% Obliczenia (boundary_type: 'dirichlet' lub 'neumann')
zad1.h = h_xt_vector(zad1.x, zad1.t, zad1.c, zad1.initial, zad1.bc_left, zad1.bc_right, 'dirichlet');

% Przygotowanie danych do animacji
zad1.points = zeros(length(zad1.t), length(zad1.x), 2);
zad1.points(:, :, 1) = repmat(zad1.x, length(zad1.t), 1);
zad1.points(:, :, 2) = zad1.h';

%% Animacja zadania 1
animate_2d(zad1.t, zad1.points, 'x', 'h(x,t)', 'Struna gitarowa');

%% Zadanie 2 - Skakanka

% Parametry zadania 2
zad2.x = 0:0.1:2;
zad2.t = 0:0.005:30;
zad2.c = 1;
zad2.L = zad2.x(end);

% Warunki brzegowe h(0, t) i h(L, t):
zad2.bc_left = @(time) sin(2 * time);
zad2.bc_right = @(time) sin(3 * time);

% Warunek początkowy h(x, 0):
zad2.initial = @(x_val) 0 * x_val;

% Obliczenia (boundary_type: 'dirichlet' lub 'neumann')
zad2.h = h_xt_vector(zad2.x, zad2.t, zad2.c, zad2.initial, zad2.bc_left, zad2.bc_right, 'dirichlet');

% Przygotowanie danych do animacji
zad2.points = zeros(length(zad2.t), length(zad2.x), 2);
zad2.points(:, :, 1) = repmat(zad2.x, length(zad2.t), 1);
zad2.points(:, :, 2) = zad2.h';

%% Animacja zadania 2
animate_2d(zad2.t, zad2.points, 'x', 'h(x,t)', 'Skakanka');

%% Zadanie 3 - Tafla jeziora

% Parametry zadania 3
zad3.x = 0:0.1:10;
zad3.t = 0:0.005:30;
zad3.c = 1;
zad3.L = zad3.x(end);

% Warunek początkowy h(x, 0):
zad3.initial = @(x_val) zad3_initial_condition(x_val, 3, 4);

% Obliczenia (boundary_type: 'dirichlet' lub 'neumann')
% Dla Neumanna bc_left/right są ignorowane, ale funkcja wymaga ich jako argumentów
zad3.h = h_xt_vector(zad3.x, zad3.t, zad3.c, zad3.initial, @(t)0, @(t)0, 'neumann');

% Przygotowanie danych do animacji
zad3.points = zeros(length(zad3.t), length(zad3.x), 2);
zad3.points(:, :, 1) = repmat(zad3.x, length(zad3.t), 1);
zad3.points(:, :, 2) = zad3.h';

%% Animacja zadania 3
animate_2d(zad3.t, zad3.points, 'x', 'h(x,t)', 'Tafla jeziora');

%% Funkcje
% Funkcja korzystająca z pętli (tylko do zadania 1.)
function [h] = zad1_h_xt_forloop(x, t, c)
    arguments
        x (1, :) {mustBeNumeric}
        t (1, :) {mustBeNumeric}
        c (1, 1) {mustBeNumeric}
    end
    L  = x(end);
    dt = t(2) - t(1);
    dx = x(2) - x(1);
    r  = c * dt / dx;
    nx = length(x);
    nt = length(t);
    h  = zeros(nx, nt);

    % Warunek początkowy t=0
    for x_idx = 1:nx
        h(x_idx, 1) = 0.5 - 0.5 * cos(2 * pi / L * x(x_idx));
    end

    % Warunek brzegowy: h(0, t) = 0 i h(L, t) = 0
    % już spełniony przez inicjalizację zerami

    % Drugi krok czasowy (t=dt)
    % Zakładamy h(x, -dt) = 0 (struna w spoczynku przed t=0)
    for x_idx = 2:nx - 1
        h(x_idx, 2) = r^2 * (h(x_idx + 1, 1) + h(x_idx - 1, 1)) + 2 * (1 - r^2) * h(x_idx, 1);
    end

    % Pozostałe kroki czasowe
    for t_idx = 2:nt - 1
        for x_idx = 2:nx - 1
            h(x_idx, t_idx + 1) = r^2 * (h(x_idx + 1, t_idx) + h(x_idx - 1, t_idx)) + 2 * (1 - r^2) * h(x_idx, t_idx) - h(x_idx, t_idx - 1);
        end
    end
end

% Funkcja uniwersalna dla równania falowego
function h = h_xt_vector(x, t, c, initial_condition, bc_left, bc_right, boundary_type)
    arguments
        x (1, :) {mustBeNumeric}
        t (1, :) {mustBeNumeric}
        c (1, 1) {mustBeNumeric}
        initial_condition {mustBeA(initial_condition, "function_handle")}
        bc_left {mustBeA(bc_left, "function_handle")}
        bc_right {mustBeA(bc_right, "function_handle")}
        boundary_type (1, :) char {mustBeMember(boundary_type, {'dirichlet', 'neumann'})} = 'dirichlet'
    end

    dt = t(2) - t(1); dx = x(2) - x(1);
    nx = length(x); nt = length(t);
    r  = c * dt / dx;
    h  = zeros(nx, nt);
    
    % Warunek początkowy (t=0)
    h(:, 1) = initial_condition(x);
    if strcmp(boundary_type, 'dirichlet')
        h([1, nx], 1) = [bc_left(t(1)), bc_right(t(1))];
    end
    
    % Drugi krok czasowy (t=dt) - Założenie: v_0 = 0 (h(t-dt) = h(t+dt))
    idx = 2:nx-1;
    h(idx, 2) = (r^2 / 2) * (h(idx+1, 1) + h(idx-1, 1)) + (1 - r^2) * h(idx, 1);
    
    if strcmp(boundary_type, 'dirichlet')
        h([1, nx], 2) = [bc_left(t(2)), bc_right(t(2))];
    else
        h(1, 2) = r^2 * h(2, 1) + (1 - r^2) * h(1, 1);
        h(nx, 2) = r^2 * h(nx-1, 1) + (1 - r^2) * h(nx, 1);
    end
    
    % Kolejne kroki czasowe (t > dt)
    for k = 2:nt-1
        h(idx, k+1) = r^2 * (h(idx+1, k) + h(idx-1, k)) + 2*(1-r^2)*h(idx, k) - h(idx, k-1);
        
        if strcmp(boundary_type, 'dirichlet')
            h([1, nx], k+1) = [bc_left(t(k+1)), bc_right(t(k+1))];
        else
            h(1, k+1) = 2*r^2*h(2, k) + 2*(1-r^2)*h(1, k) - h(1, k-1);
            h(nx, k+1) = 2*r^2*h(nx-1, k) + 2*(1-r^2)*h(nx, k) - h(nx, k-1);
        end
    end
end

% Warunek początkowy dla zadania 3
function h = zad3_initial_condition(x_val, a, b)
    h = zeros(size(x_val));
    
    mask = (x_val >= a) & (x_val <= b);
    n_disturb = sum(mask);
    
    if n_disturb > 0
        h(mask) = -0.5 + 0.5 * cos(linspace(0, 2*pi, n_disturb));
    end
end

% Funkcja animacji
function animate_2d(t, points, x_label, y_label, title_str)
arguments
    t         (1, :)    {mustBeNumeric}
    points    (:, :, 2) {mustBeNumeric}
    x_label   (1, :)    char = 'x'
    y_label   (1, :)    char = 'y'
    title_str (1, :)    char = ''
end

    dt = t(2) - t(1);

    figure;
    hLine = plot(squeeze(points(1, :, 1)), squeeze(points(1, :, 2)), 'b-', 'LineWidth', 2);
    xlim([min(points(:, :, 1), [], 'all'), max(points(:, :, 1), [], 'all')]);
    ylim([min(points(:, :, 2), [], 'all')-0.1, max(points(:, :, 2), [], 'all')+0.1]);
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
