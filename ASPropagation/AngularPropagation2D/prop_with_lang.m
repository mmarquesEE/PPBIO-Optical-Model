clear all
close all
clc

t0 = cputime;

%% Physical Constants and Parameters
lambda0_um = 0.67;
lambda0 = lambda0_um*1e-6;
%% Material properties
% Load refractive index data
gold_refr_idx_table = readtable("gold.txt");
lambda_gold = gold_refr_idx_table.lambda;
n_gold = gold_refr_idx_table.n;
k_gold = gold_refr_idx_table.k;

water_refr_idx_table = readtable("water.txt");
lambda_water = water_refr_idx_table.lambda;
n_water = water_refr_idx_table.n;
k_water = water_refr_idx_table.k;

% Calculate refractive indices
% N-BK7
n1 = sqrt(...
    1+1.03961212./(1 - 0.00600069867./lambda0_um.^2) + ...
    0.231792344./(1 - 0.0200179144./lambda0_um.^2) + ...
    1.01046945./(1 - 103.560653./lambda0_um.^2) ...
);

% Gold
n2 = spline(lambda_gold, n_gold, lambda0_um) + ...
    1j*spline(lambda_gold, k_gold, lambda0_um);

% Water
n3 = spline(lambda_water, n_water, lambda0_um) + ...
    1j*spline(lambda_water, k_water, lambda0_um);

% thickness of gold layer
d_tf = 50e-9;

lambda1 = lambda0/n1;
%% Gaussian Beam Parameters
theta01 = 10/180*pi; % Divergence angle [rad]
theta02 = 0.1/180*pi;

k = 2*pi/lambda1; % Wavevector

% Beam 1 parameters
w01 = lambda1/(pi*theta01);
z01 = lambda1/(pi*theta01^2);

% Beam 2 parameters
w02 = lambda1/(pi*theta02);
z02 = lambda1/(pi*theta02^2);
%% Angular Spectrum Propagation
d_alpha = 0.02/180*pi;
ds = d_alpha;

dz = d_tf/4;

s = ds*linspace(-0.5,0.5,2^13);
%% Incident Angle Configuration
theta = 68.8/180*pi;

% Coordinate transformation
x = s*cos(theta);
z = -s*sin(theta);
%% Initialize Fields
U1 = gb_complex_amplitude(0.3, x, z, k, w01, z01);
U2 = gb_complex_amplitude(0.3, x, z, k, w02, z02);
% Precompute angular spectra
[Ai1, fs1] = angular_spectrum(U1, s);
[Ai2, fs2] = angular_spectrum(U2, s);
%% Using R-TMM
% Langmuir adsorption parameters (Biotin example)
k_on = 160;       % [M⁻¹s⁻¹]
k_off = 0.054;     % [s⁻¹]
C = 0.82e-3;         % [M]

% Simulation time parameters
t_total = 40;      % [s] Total simulation time
num_frames = 10;    % Number of animation frames

%% Solve Langmuir Kinetics
[time_points, theta_frac] = ode45(@(t,y) k_on*C*(1.2-y) - k_off*y, ...
                                 linspace(0, t_total, num_frames), 0);
diff_n = 1e-3*theta_frac'; % Refractive index change

%% Initialize Figures
fig1 = figure; % Bottom view
fig2 = figure('Position', [100, 100, 1200, 150]);% Field Distribution
fig3 = figure('Position', [300, 300, 800, 400]); % Sensogram
fig4 = figure; % Far Field
ax1 = axes(fig1);
title(ax1, 'Beam Cross-Section');
ax2 = axes(fig2);
colormap(ax2, parula(256));
ax3 = axes(fig3);
cmap = parula(num_frames); % Colormap for time evolution
ax4 = axes(fig4);

h_sensogram = plot(ax3, NaN, NaN, 'b-', 'LineWidth', 2);

zt = ((-7e-6):dz:1e-6)';

filename1 = 'bottom_view_biotin.gif';
filename2 = 'side_view_biotin.gif';
filename3 = 'Sensorgram_biotin.gif';
filename4 = 'FarField_biotin.gif';

delay_time = 0.1;
for t = 1:length(diff_n)
    % Calculate transmitted field with current refractive index
    [At1, rp1, ~] = tmm_3p_TM(Ai1, fs1, zt, [n1;n2;n3 + diff_n(t)], d_tf, lambda0);
    [At2, rp2, ~] = tmm_3p_TM(Ai2, fs1, zt(1:2), [n1;n2;1], d_tf, lambda0);
    % Reconstruct fields
    [Ut1, sp1] = i_angular_spectrum(At1, fs1);
    [Ut2, ~] = i_angular_spectrum(At2, fs2);
    
    [Ux, Uy] = meshgrid(Ut1(1,:), Ut2(1,:));
    [XX , ] = meshgrid(x,x);
    
    idx = (XX > 0.4e-5) & (XX < 1.2e-5);
    %% Plot 1: Beam Cross-Section
    imagesc(ax1,sqrt(abs(Ux(idx)).*abs(Uy(idx)))')
    title(ax1, sprintf('t = %.1f s (Δn = %.2e)', time_points(t), diff_n(t)));
    colormap(ax1, 'gray'); % Set to grayscale
    xlabel off
    ylabel off
    %% Plot 2: Field Distribution
    hold on
    imagesc(ax2,-1e6*sp1, 1e6*zt, abs(Ut1))
    title(ax2, ['$|E_{xz}|$ ($\lambda=',num2str(lambda0_um,2),...
    '\mu$m, $\theta = ',num2str(theta/pi*180,2),'$deg, ',...
    '$\theta_{||}=',num2str(theta01/pi*180),'$deg)'],...
    Interpreter='latex', FontSize=16)
    xlabel('x($\mu$m)', Interpreter='latex', FontSize=16)
    ylabel('z($\mu$m)', Interpreter='latex', FontSize=16);
    xlim(ax2,40*[-1, 1])
    set(ax2, 'YDir', 'normal')  % <<< This line inverts the Z axis
    hold off
    %% Update Sensogram
    set(h_sensogram, 'XData', time_points(1:t), 'YData', diff_n(1:t));
    xlabel(ax3, 'Time (s)');
    ylabel(ax3, 'RI');
    title(ax3, 'Sensogram: RI over Time');

   % Compute reflected angular spectra
    Ar1 = rp1 .* Ai1; % Reflected spectrum for beam 1
    
    % Apply propagation to Z_cam using angular spectrum method
    k = 2*pi / lambda1; % Wavevector in medium n1
    
    % Propagation phase factor for each spatial frequency component
    Z_cam  = -5e-6; % Distance to camera plane [meters]
    H1 = exp(1j * k * Z_cam * sqrt(1 - (lambda1 * fs1).^2));
    prop_mask = real(sqrt(1 - (lambda1 * fs1).^2)) > 0;
    H1 = H1 .* prop_mask;
    Ar1_propagated = Ar1 .* H1;
    
    % Compute fields at Z_cam
    [Ufar1, x_cam1] = i_angular_spectrum(Ar1_propagated, fs1);
    
    % Calculate the intensity (product of intensities from both beams)
    Intensity = sqrt(abs(Ufar1).*abs(Ufar1'));
    Intensity = Intensity / max(Intensity(:));

    %% Plot Orthogonal Plane to Reflected Beam (Camera-like)

    % === Step 1: Reflected beam direction vector ===
    theta_r = theta;  % Reflected angle equals incident angle in planar geometry
    n_vec = [sin(theta_r); 0; -cos(theta_r)]; % Reflected beam direction
    
    % === Step 2: Vectors in the orthogonal plane ===

    v1 = [0; 1; 0]; v1 = v1 / norm(v1); % Tangent vector 1
    v2 = cross(n_vec, v1);       v2 = v2 / norm(v2); % Tangent vector 2
    
     % Create grid with proper coordinate ranges
    Lx = 10e-6; % X-span (along v2 direction)
    Ly = 100e-6; % Y-span (along v1 direction)
    N = 600;
    
    % X: 0 to positive, Y: symmetric negative/positive
    [u, v] = meshgrid(linspace(0, Lx, N), linspace(-Ly, Ly, N));
    origin = [0*1e-6,0*1e-6,-10*1e-6];%n_vec * Z_cam;  % Center of the plane at Z_cam
    
    X_plane = origin(1) + u*v2(1) + v*v1(1);
    Y_plane = origin(2) + u*v2(2) + v*v1(2);
    Z_plane = origin(3) + u*v2(3) + v*v1(3);

    % === Step 4: Interpolate intensity onto rotated plane ===
    % Use existing camera grid (x_cam1, x_cam2), Intensity
    [cam_X, cam_Y] = meshgrid(x_cam1, x_cam1);  % Original far-field grid
    
    % Interpolate with safety (fill missing with 0)
    RotatedIntensity = interp2(cam_X, cam_Y, Intensity, X_plane, Y_plane, 'linear', 0);
    % === Step 5: Plot the rotated camera plane ===
    surf(ax4, X_plane*1e6, Y_plane*1e6, Z_plane*1e6, RotatedIntensity, ...
     'EdgeColor', 'none', 'FaceAlpha', 1);
    xlabel('X (um)');ylabel('Y (um)');zlabel('Z (um)');
    title(ax4, sprintf('Orthogonal Plane View (t = %.1f s)', time_points(t)));
    colormap(ax4,"gray");
    colorbar;
    % Beam direction vector overlay
    hold on
    quiver3(ax4, 0, 0, 0, n_vec(1), n_vec(2), n_vec(3), 1e-5, ...
            'Color', 'cyan', 'LineWidth', 2, 'MaxHeadSize', 2);
    text(n_vec(1)*1.1e-5, n_vec(2)*1.1e-5, n_vec(3)*1.1e-5, 'Beam →', ...
         'Color', 'cyan', 'FontSize', 12);

    % Captura o frame e converte para imagem indexada
    frame1 = getframe(fig1);
    im1 = frame2im(frame1);
    [imind1, cm1] = rgb2ind(im1, 256);

    % Captura o frame e converte para imagem indexada
    frame2 = getframe(fig2);
    im2 = frame2im(frame2);
    [imind2, cm2] = rgb2ind(im2, 256);

    % Captura o frame e converte para imagem indexada
    frame3 = getframe(fig3);
    im3 = frame2im(frame3);
    [imind3, cm3] = rgb2ind(im3, 256);

    % Captura o frame e converte para imagem indexada
    frame4 = getframe(fig4);
    im4 = frame2im(frame4);
    [imind4, cm4] = rgb2ind(im4, 256);
    
    % Escreve no GIF
    if t == 1
        imwrite(imind1, cm1, filename1, 'gif', 'Loopcount', inf, 'DelayTime', delay_time);
        imwrite(imind2, cm2, filename2, 'gif', 'Loopcount', inf, 'DelayTime', delay_time);
        imwrite(imind3, cm3, filename3, 'gif', 'Loopcount', inf, 'DelayTime', delay_time);
        imwrite(imind4, cm4, filename4, 'gif', 'Loopcount', inf, 'DelayTime', 0.1);

    else
        imwrite(imind1, cm1, filename1, 'gif', 'WriteMode', 'append', 'DelayTime', delay_time);
        imwrite(imind2, cm2, filename2, 'gif', 'WriteMode', 'append', 'DelayTime', delay_time);
        imwrite(imind3, cm3, filename3, 'gif', 'WriteMode', 'append', 'DelayTime', delay_time);
        imwrite(imind4, cm4, filename4, 'gif', 'WriteMode', 'append', 'DelayTime', delay_time);

    end
end
%% Plot
execution_time = cputime - t0



function [E, rp, tp] = tmm_3p_TM(E_TM, fx, z, n, d, lambda0)
% Transfer Matrix Method for a 3phase system, TM mode

lambda = lambda0./n;
arg = 1./lambda(1).^2 - fx.^2;

fz1 = sqrt(arg);
fz1(arg < 0) = 1j*fz1(arg < 0);
fz = (n./n(1)).*fz1;

theta = real(asin(fx*lambda(1)));

xi = sqrt(n.^2 - n(1)^2*sin(theta).^2);
q = xi./n.^2;
beta = 2*pi*d/lambda0*xi;

m11= cos(beta(2));
m12 = -1j/q(2)*sin(beta(2));
m21 = -1j*q(2)*sin(beta(2));
m22 = cos(beta(2));

rp = ((m11 + m12.*q(end,:)).*q(1,:) - (m21 + m22.*q(end,:)))./...
     ((m11 + m12.*q(end,:)).*q(1,:) + (m21 + m22.*q(end,:)));

tp = 2*q(1)./((m11 + m12.*q(end)).*q(1) + (m21 + m22.*q(end)));

z1 = z(z <= 0);
z2 = z((z > 0) & (z <= d));
z3 = z(z > d) - d;


p1 = 2*pi*z1*fz(1,:);
p2 = 2*pi*z2*fz(2,:);
p3 = 2*pi*z3*fz(3,:);

Ht = n(1)*E_TM;

E1 = -q(1)*Ht.*(exp(1j*p1) - rp.*exp(-1j*p1));
E2 = -q(2)*Ht.*(exp(1j*p2) - rp.*exp(-1j*p2));
E3 = -q(end)*tp.*Ht.*exp(1j*p3);

E = [E1;E2;E3];

end
