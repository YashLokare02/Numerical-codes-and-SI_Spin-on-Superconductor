clc;
clear;
close all;

load('Delta.mat');

T = [1.08 1.25 1.5 1.75 2 2.25 2.5 2.75 3 3.25 3.5 3.75 4]; 

figure(1)
hold on
plot(T, Gamma_mv, 'b');
title('Plot for the thermal broadening parameter vs. temperature', 'FontSize', 20);
xlabel('Temperature [K]', 'FontSize', 18);
ylabel('$\Gamma$', 'Interpreter', 'latex', 'FontSize', 18);
hold off

figure(2)
hold on
plot(T, delta_mv, 'b');
title('Plot for the superconducting energy gap vs. temperature', 'FontSize', 20);
xlabel('Temperature [K]', 'FontSize', 18);
ylabel('$\Delta$', 'Interpreter', 'latex', 'FontSize', 18);
hold off