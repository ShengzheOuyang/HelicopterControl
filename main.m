clc
clear
% calculate Arbeitpunkt
run("cal_xs.m")  

% Linearisierung
run('solve.m')

% calculate L for Luenberger Beobachter
run("Beobachter.m")
