%{
Systemy Cyber-Fizyczne - Laboratorium 1
Elementy mechaniki klasycznej w modelowaniu układów mechanicznych
%}

close all
clear
clc

% Zadanie #1

%{
Wyprowadzenie wzoru na transmitancję:
    ma = sum(F)
    m*ddot_y = -b(dot_y - dot_u) - k(y - u)
    m*ddot_y + b*dot_y + k*y = b*dot_u + k*u
    ->LAPLACE
    Y(s) * (m*s^2 + b*s + k) = U(s) * (b*s + k)
    G(s) = Y(s)/U(s) = (b*s + k)/(m*s^2 + b*s + k)
%}

z1.t = 0:0.01:150;
z1.k = 0.1; z1.b = 0.1; z1.m = 1;
z1.sys = tf([z1.b z1.k], [z1.m z1.b z1.k]);
z1.sys

figure;
set(gcf, 'Position', [100, 100, 1200, 350]);

t = tiledlayout(1, 3, 'TileSpacing','compact','Padding','compact');
title(t, 'Zadanie 1.')

nexttile
plot(z1.t, step(z1.sys, z1.t, RespConfig("Amplitude",1)), "LineWidth",1)
title("Skok A=1")
xlabel("Czas t"), ylabel("y(t)")
grid on

nexttile
plot(z1.t, impulse(z1.sys, z1.t, RespConfig("Amplitude",1)), "LineWidth",1)
title("Impuls A=1")
xlabel("Czas t")
grid on

nexttile
plot(z1.t, impulse(z1.sys, z1.t, RespConfig("Amplitude",5)), "LineWidth",1)
title("Impuls A=5")
xlabel("Czas t")
grid on

% Zadanie #2

%{
Wyprowadzenie wzoru na transmitancję:
    ma = sum(F)
    m*ddot_y = -b*dot_y - k*y + u
    m*ddot_y + b*dot_y + k*y = u
    ->LAPLACE
    Y(s) * (m*s^2 + b*s + k) = U(s)
    G(s) = Y(s)/U(s) = 1/(m*s^2 + b*s + k)
%}

z2.t = 0:0.01:150;
z2.k = 0.1; z2.b = 0.1; z2.m = 1;
z2.sys = tf(1, [z2.m z2.b z2.k]);
z2.sys

figure;
set(gcf, 'Position', [100, 100, 1200, 350]);

t = tiledlayout(1,2, 'TileSpacing','compact','Padding','compact');
title(t, 'Zadanie 2.')

nexttile
plot(z2.t, step(z2.sys, z2.t, RespConfig("Amplitude",5)), "LineWidth",1)
title("Skok A=5")
xlabel("Czas t"), ylabel("y(t)")
grid on

nexttile
plot(z2.t, impulse(z2.sys, z2.t, RespConfig("Amplitude",3)), "LineWidth",1)
title("Impuls A=3")
xlabel("Czas t")
grid on

% Zadanie #3

%{
Wyprowadzenie wzoru na transmitancję:
    ma = sum(F)
    {
        m1*ddot_y1 = -k1(y1 - u) + k2(y2 - y1) + b(dot_y2 - dot_y1)
        m2*ddot_y2 = -k2(y2 - y1) - b(dot_y2 - dot_y1)
    }
    {
        Y1(s) * (m1*s^2 + b*s + (k1 + k2)) + Y2(s) * (-b*s - k2) = U(s) * k1
        Y2(s) * (m2*s^2 + b*s + k2) + Y1(s) * (-b*s - k2) = 0
    }
%}

syms s Y1 Y2 U m1 m2 k1 k2 b;
z3.sys_sym = simplify(solve( ...
        [Y1 * (m1*s^2 + b*s + (k1 + k2)) + Y2 * (-b*s - k2) == U * k1, Y2 * (m2*s^2 + b*s + k2) + Y1 * (-b*s - k2) == 0], ...
        [Y1, Y2]).Y2/U);
[z3.sys_sym_num, z3.sys_sym_den] = numden(z3.sys_sym);
z3.sys_sym_num = collect(expand(z3.sys_sym_num), s);
z3.sys_sym_den = collect(expand(z3.sys_sym_den), s);
z3.sys_sym = z3.sys_sym_num/z3.sys_sym_den;
disp(z3.sys_sym);

z3.t = 0:0.001:5;
z3.m1 = 150; z3.m2 = 1500; z3.k1 = 200000; z3.k2 = 48000; z3.b = 10000;

figure;
set(gcf, 'Position', [100, 100, 1200, 500]);

t = tiledlayout(2,2, 'TileSpacing','compact','Padding','compact');
title(t, 'Zadanie 3.')

z3.sys_num = subs(coeffs(z3.sys_sym_num, s), [m1,m2,k1,k2,b], [z3.m1,z3.m2,z3.k1,z3.k2,z3.b]);
z3.sys_den = subs(coeffs(z3.sys_sym_den, s), [m1,m2,k1,k2,b], [z3.m1,z3.m2,z3.k1,z3.k2,z3.b]);
z3.sys_num = double(z3.sys_num(end:-1:1));
z3.sys_den = double(z3.sys_den(end:-1:1));

z3.sys = tf(z3.sys_num, z3.sys_den);

% --- Skok A=4
nexttile
plot(z3.t, step(z3.sys, z3.t, RespConfig("Amplitude",4)), "LineWidth",1)
title("Skok A=4")
xlabel("Czas t"), ylabel("y(t)")
grid on

% --- Skok A=-2
nexttile
plot(z3.t, step(z3.sys, z3.t, RespConfig("Amplitude",-2)), "LineWidth",1)
title("Skok A=-2")
xlabel("Czas t")
grid on

b_vals = [0, 50, 100, 500, 1000, 10000, 50000, 100000];
colors = turbo(length(b_vals));

z3.t = 0:0.001:3;
% --- Impuls vs b
nexttile
hold on
for i = 1:length(b_vals)
    z3.b = b_vals(i);

    z3.sys_num = subs(coeffs(z3.sys_sym_num, s), [m1,m2,k1,k2,b], [z3.m1,z3.m2,z3.k1,z3.k2,z3.b]);
    z3.sys_den = subs(coeffs(z3.sys_sym_den, s), [m1,m2,k1,k2,b], [z3.m1,z3.m2,z3.k1,z3.k2,z3.b]);
    z3.sys_num = double(z3.sys_num(end:-1:1));
    z3.sys_den = double(z3.sys_den(end:-1:1));

    z3.sys = tf(z3.sys_num, z3.sys_den);

    plot(z3.t, impulse(z3.sys, z3.t, RespConfig("Amplitude",4)), "LineWidth",1, ...
        'Color',colors(i,:), 'DisplayName',sprintf("b=%d", z3.b))
end
title("Impuls A=4 - wpływ b")
xlabel("Czas t")
grid on
legend('Location','southoutside', 'Orientation','horizontal', 'NumColumns',4)

z3.b = 10000;
k2_vals = [0, 50, 100, 500, 1000, 10000, 50000, 100000];
colors = turbo(length(k2_vals));

% --- Impuls vs k2
nexttile
hold on
for i = 1:length(k2_vals)
    z3.k2 = k2_vals(i);

    z3.sys_num = subs(coeffs(z3.sys_sym_num, s), [m1,m2,k1,k2,b], [z3.m1,z3.m2,z3.k1,z3.k2,z3.b]);
    z3.sys_den = subs(coeffs(z3.sys_sym_den, s), [m1,m2,k1,k2,b], [z3.m1,z3.m2,z3.k1,z3.k2,z3.b]);
    z3.sys_num = double(z3.sys_num(end:-1:1));
    z3.sys_den = double(z3.sys_den(end:-1:1));

    z3.sys = tf(z3.sys_num, z3.sys_den);

    plot(z3.t, impulse(z3.sys, z3.t, RespConfig("Amplitude",4)), "LineWidth",1, ...
        'Color',colors(i,:), 'DisplayName',sprintf("k_2=%d",z3.k2))
end
title("Impuls A=4 - wpływ k_2")
xlabel("Czas t")
grid on
legend('Location','southoutside', 'Orientation','horizontal', 'NumColumns',4)