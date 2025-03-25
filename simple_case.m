clear all
clc

% Inclua o pacote de simbolos
syms lambda theta0 z f_c f_f n_p

W0 = lambda/(pi*theta0);
z0 = lambda/(pi*theta0^2);

z0c = (f_c^2/(2*z0^2))*z0;
z0f = z0c/(1 + (z0c/f_f)^2);
z0r = n_p*z0f;

