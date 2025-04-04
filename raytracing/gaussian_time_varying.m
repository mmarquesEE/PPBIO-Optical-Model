clc; clear; close all;

%% Parâmetros do feixe
lambda = 0.4e-6;  % Comprimento de onda (m)
w0 = 0.1e-6;     % Raio de cintura do feixe (m)
z0 = pi * w0^2 / lambda;  % Distância de Rayleigh (m)
k = 2 * pi / lambda;  % Número de onda
c = 3e8;  % Velocidade da luz no vácuo (m/s)
nu = c / lambda;  % Frequência temporal
omega = 2 * pi * nu;  % Frequência angular
A0 = 1;  % Normalização da amplitude
fps = 30;  % Frames por segundo
num_frames = 60;  % Duração da animação
delay_time = 0.1; % Tempo entre frames no GIF (segundos)

%% Funções do feixe Gaussiano
w = @(z) w0 * sqrt(1 + (z ./ z0).^2);
R = @(z) z .* (1 + (z0 ./ z).^2);
psi = @(z) atan(z ./ z0);

%% -------- Gráfico (a) - Cintura do Feixe --------
% z_vals = linspace(-2*z0, 2*z0, 200);
% w_vals = w(z_vals);
% 
% figure;
% hold on;
% fill([z_vals, fliplr(z_vals)], [-w_vals, fliplr(w_vals)], 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
% plot(z_vals, w_vals, 'r', 'LineWidth', 1.5);
% plot(z_vals, -w_vals, 'r', 'LineWidth', 1.5);
% plot([0, 0], [-3*w0, 3*w0], 'k--', 'LineWidth', 0.8); % Linha central z = 0
% plot([-2*z0, 2*z0], [0, 0], 'k--', 'LineWidth', 0.8); % Eixo x
% 
% xlabel('z (m)');
% ylabel('x (m)');
% title('Perfil da Cintura do Feixe Gaussiano');
% xlim([-2*z0, 2*z0]);
% ylim([-3*w0, 3*w0]);
% grid on;
% hold off;

%% -------- Criando o GIF do Campo Elétrico e Frentes de Onda --------
z_vals = linspace(-2 * z0, 2 * z0, 50);
x_vals = linspace(-3 * w0, 3 * w0, 50);
[X, Z] = meshgrid(x_vals, z_vals);
rho2 = X.^2;  % y = 0 para simplificação

filename = 'electric_field.gif';

for t = linspace(0, 2*pi/omega, num_frames)
    clf;
    hold on;
    
    % Amplitude complexa com dependência temporal
    U = A0 * (w0 ./ w(Z)) .* exp(-rho2 ./ w(Z).^2) ...
        .* exp(-1j * (k * Z + k * rho2 ./ (2 * R(Z)) - psi(Z) - omega * t));

    % Vetor de onda
    kvx =  (k*X)./(Z.*(z0^2./Z.^2 + 1));
    kvz = k - 1./(z0.*(Z.^2./z0.^2 + 1)) - ...
        (k.*X.^2)./(2.*Z.^2.*(z0.^2./Z.^2 + 1)) + ...
        (k.*X.^2.*z0.^2)./(Z.^4.*(z0.^2./Z.^2 + 1).^2);

    % Campo elétrico vetorial
    E0 = 1;  
    Ex = -E0 * U;  
    Ez = E0 * (X ./ (Z + 1j * z0)) .* U;
    kvx = kvx.*sqrt(real(Ex).^2 + real(Ez).^2);
    kvz = kvz.*sqrt(real(Ex).^2 + real(Ez).^2);

    % Máscara para limitar os vetores à cintura do feixe
    mask = (abs(X) <= w(Z));  
    Ex(~mask) = NaN;
    Ez(~mask) = NaN;
    kvx(~mask) = NaN;
    kvz(~mask) = NaN;

    quiver(Z, X, real(Ez), real(Ex), 'r', 'AutoScale', 'on', 'AutoScaleFactor', 1.5);
    quiver(Z, X, kvz, kvx, 'b', 'AutoScale', 'on', 'AutoScaleFactor', 1.5);

    plot([0, 0], [-3*w0, 3*w0], 'k--', 'LineWidth', 0.8); % Linha central
    plot([-2*z0, 2*z0], [0, 0], 'k--', 'LineWidth', 0.8); % Eixo x
%     plot([-2*z0, 2*z0], [-3*w0, 3*w0], 'k--', 'LineWidth', 0.8); 

    axis equal

    xlabel('z (m)');
    ylabel('x (m)');
    title('Oscilação do Campo Elétrico e Frentes de Onda');
    xlim([-2*z0, 2*z0]);
    ylim([-3*w0, 3*w0]);
    
    grid on;
    hold off;
    
    % Captura o frame e converte para imagem indexada
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    
    % Escreve no GIF
    if t == 0
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', delay_time);
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', delay_time);
    end
end