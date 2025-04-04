close all
clear all
clc

%% Parâmetros do sistema



lambda = 670e-9;
theta0x = 5*0.5*pi/180;
theta0y = 30*0.5*pi/180;

D = 1e-2;

fc = (D/2)/tan(theta0y);
ff = 20e-3;


