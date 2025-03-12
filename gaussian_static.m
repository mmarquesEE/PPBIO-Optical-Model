clc; clear; close all;

%% Parâmetros do feixe
lambda = 0.7e-6;  % Comprimento de onda (m)
w0 = 0.37e-6;        % Raio de cintura do feixe (m)
z0 = pi*w0^2/lambda;        % Distância de Rayleigh (m)
k = 2 * pi / lambda;  % Número de onda
A0 = 1;  % Normalização da amplitude

%% Funções do feixe Gaussiano
w = @(z) w0 * sqrt(1 + (z ./ z0).^2);
R = @(z) z .* (1 + (z0 ./ z).^2);
psi = @(z) atan(z ./ z0);

%% -------- Gráfico (a) - Cintura do Feixe --------
z_vals = linspace(-2*z0, 2*z0, 200);
w_vals = w(z_vals);

figure;
hold on;
fill([z_vals, fliplr(z_vals)], [-w_vals, fliplr(w_vals)], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
plot(z_vals, w_vals, 'r', 'LineWidth', 1.5);
plot(z_vals, -w_vals, 'r', 'LineWidth', 1.5);
plot([0, 0], [-3*w0, 3*w0], 'k--', 'LineWidth', 0.8); % Linha central z = 0
plot([-2*z0, 2*z0], [0, 0], 'k--', 'LineWidth', 0.8); % Eixo x

xlabel('z (m)');
ylabel('x (m)');
title('Perfil da Cintura do Feixe Gaussiano');
xlim([-2*z0, 2*z0]);
ylim([-3*w0, 3*w0]);
grid on;
hold off;

%% -------- Gráfico (b) - Campo Elétrico e Frentes de Onda --------
z_vals = linspace(-2 * z0, 2 * z0, 50);
x_vals = linspace(-3 * w0, 3 * w0, 30);  % Maior resolução na área central
[X, Z] = meshgrid(x_vals, z_vals);
rho2 = X.^2;  % y = 0 para simplificação

% Amplitude complexa U(x,z)
U = A0 * (w0 ./ w(Z)) .* exp(-rho2 ./ w(Z).^2) ...
    .* exp(-1j * (k * Z + k * rho2 ./ (2 * R(Z)) - psi(Z)));

% Campo elétrico vetorial
E0 = 1;  
Ex = -E0 * U;  
Ez = E0 * (X ./ (Z + 1j * z0)) .* U; 

% Máscara para limitar os vetores à cintura do feixe
mask = abs(X) <= w(Z);  
Ex(~mask) = NaN;  % Remove valores fora da cintura do feixe
Ez(~mask) = NaN;


figure;
hold on;
quiver(Z, X, real(Ez), real(Ex), ...
    'r', 'AutoScale', 'on', 'AutoScaleFactor', 1);
plot([0, 0], [-3*w0, 3*w0], 'k--', 'LineWidth', 0.8); % Linha central
plot([-2*z0, 2*z0], [0, 0], 'k--', 'LineWidth', 0.8); % Eixo x

xlabel('z (m)');
ylabel('x (m)');
title('Vetores do Campo Elétrico e Frentes de Onda');
xlim([-2*z0, 2*z0]);
ylim([-3*w0, 3*w0]);
grid on;
hold off;
