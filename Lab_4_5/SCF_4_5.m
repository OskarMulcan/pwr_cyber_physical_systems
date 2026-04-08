%{
Systemy Cyber-Fizyczne - Laboratorium 4 i 5
Analiza modeli fizycznych w kontekście własności w przestrzeni fazowej
Przestrzeń fazowa a właściwości obiektów - wybrane modele obiektów
%}

close all; clc; clear;

% Zadanie #1
zad1.dot_x = @(t, x) 3 * exp(-t);
zad1.t = 0:0.1:5;
zad1.x_0 = 0;
[zad1.t, zad1.x] = ode45(zad1.dot_x, zad1.t, zad1.x_0);

figure
plot(zad1.t, zad1.x);

% Zadanie #2
zad2.dot_x = @(t, x) [x(2); -x(2) + 4*x(1) + sin(10*t)];
zad2.t = 0:0.1:5;
zad2.x_0 = [0; 0];
[zad2.t, zad2.x] = ode45(zad2.dot_x, zad2.t, zad2.x_0);

figure
plot(zad2.t, zad2.x(:,1));

% Zadanie #3
function zad3_phase_portrait(A, T, quiver_range_x1, quiver_range_x2, traj_range_x1, traj_range_x2)
    arguments
        A               (2,2) {mustBeNumeric}
        T               (1,1) {mustBeNumeric, mustBePositive}
        quiver_range_x1 (1,:) {mustBeNumeric}
        quiver_range_x2 (1,:) {mustBeNumeric}                  = quiver_range_x1
        traj_range_x1   (1,:) {mustBeNumeric}                  = quiver_range_x1
        traj_range_x2   (1,:) {mustBeNumeric}                  = quiver_range_x2
    end
    dot_x = @(t, x) A * x;
    t = 0:0.01:T;

    [tmp1, tmp2] = meshgrid(quiver_range_x1, quiver_range_x2);
    x_0 = [tmp1(:), tmp2(:)];
    dot_x_0 = (dot_x(t,x_0'))';

    figure
    quiver(x_0(:,1), x_0(:,2), dot_x_0(:,1), dot_x_0(:,2))
    hold on

    [tmp1, tmp2] = meshgrid(traj_range_x1, traj_range_x2);
    for traj_x_0 = [tmp1(:), tmp2(:)]'
        [~, traj] = ode45(dot_x, t, traj_x_0);
        plot(traj(:,1), traj(:,2));
    end

    xlim([min(x_0(:,1)) max(x_0(:,1))]);
    ylim([min(x_0(:,2)) max(x_0(:,2))]);
    xlabel('x_1'); ylabel('x_2');
    title(sprintf(['A=', mat2str(A)]));
    grid on;
end

%a)
A = [-1 4;
     -2 5];
zad3_phase_portrait(A, 1, -2.8:0.4:2.8)

%b)
A = [-3  0;
      0 -2];
zad3_phase_portrait(A, 2, -8:2:8)

%c)
A = [4  0;
     2 -1];
zad3_phase_portrait(A, 1, -2.8:0.4:2.8)

%d)
A = [2   -3;
     1/3  4];
zad3_phase_portrait(A, 1, -2.8:0.4:2.8)

%e)
A = [-7  1;
     -4 -3];
zad3_phase_portrait(A, 2, -8:2:8)

%f)
A = [-2  3;
     -3 -2];
zad3_phase_portrait(A, 2, -8:2:8)

%g)
A = [ 2  3;
     -3  2];
zad3_phase_portrait(A, 1, -2.8:0.4:2.8)

%h)
A = [ 0  1;
     -5  0];
zad3_phase_portrait(A, 2, -8:2:8)

