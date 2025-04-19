clear all
close all
clc

t0 = cputime;

%% Vaccum Wavelength

lambda0_um = 0.67;
lambda0 = lambda0_um*1e-6;

%% Materials

gold_refr_idx_table = readtable("gold.txt");
lambda_gold = gold_refr_idx_table.lambda;
n_gold = gold_refr_idx_table.n;
k_gold = gold_refr_idx_table.k;

water_refr_idx_table = readtable("water.txt");
lambda_water = water_refr_idx_table.lambda;
n_water = water_refr_idx_table.n;
k_water = water_refr_idx_table.k;

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

%% Laser (defined at the substrate)

theta01 = 10/180*pi;
theta02 = 0.1/180*pi;

k = 2*pi/lambda1;

w01 = lambda1/(pi*theta01);
z01 = lambda1/(pi*theta01^2);

w02 = lambda1/(pi*theta02);
z02 = lambda1/(pi*theta02^2);

%% Angular Spectrum Propagation

d_alpha = 0.02/180*pi;
ds = d_alpha;

dz = d_tf/4;

s = ds*linspace(-0.5,0.5,2^13);

%%
theta = 68.8/180*pi;

x = s*cos(theta);
z = -s*sin(theta);

U1 = gb_complex_amplitude(0.3, x, z, k, w01, z01);
U2 = gb_complex_amplitude(0.3, x, z, k, w02, z02);

[Ai1, fs1] = angular_spectrum(U1, s);
[Ai2, fs2] = angular_spectrum(U2, s);

%% Using R-TMM

diff_n = linspace(-1e-2,3e-2,60);

fig1 = figure;
fig2 = figure('Position', [100, 100, 1200, 150]);

ax1 = axes(fig1);
ax2 = axes(fig2);

zt = ((-7e-6):dz:1e-6)';

filename1 = 'bottom_view.gif';
filename2 = 'side_view.gif';
delay_time = 0.1;

for t = 1:length(diff_n)
    t
    At1 = tmm_3p_TM(Ai1, fs1, zt, [n1;n2;n3 + diff_n(t)], d_tf, lambda0);
    At2 = tmm_3p_TM(Ai2, fs1, zt(1:2), [n1;n2;1], d_tf, lambda0);
    
    [Ut1, sp1] = i_angular_spectrum(At1, fs1);
    [Ut2, sp2] = i_angular_spectrum(At2, fs2);
    
    [Ux, Uy] = meshgrid(Ut1(1,:), Ut2(1,:));
    [XX , ] = meshgrid(x,x);
    
    idx = (XX > 0.4e-5) & (XX < 1.2e-5);
    %%
    imagesc(ax1,sqrt(abs(Ux(idx)).*abs(Uy(idx)))')
    xlabel off
    ylabel off

    %%

    hold on
    imagesc(ax2,-1e6*sp1, 1e6*zt, abs(Ut1))
%     plot(1e6*[min(sp1),max(sp1)], [0,0], 'r', 1e6*[min(sp1),max(sp1)], 1e6*d_tf*[1,1], 'r');
    hold off
    
    title(['$|E_{xz}|$ ($\lambda=',num2str(lambda0_um,2),...
        '\mu$m, $\theta = ',num2str(theta/pi*180,2),'$deg, ',...
        '$\theta_{||}=',num2str(theta01/pi*180),'$deg)'],...
        Interpreter='latex', FontSize=16)
    xlabel('x($\mu$m)', Interpreter='latex', FontSize=16)
    ylabel('z($\mu$m)', Interpreter='latex', FontSize=16);
%     axis equal
    xlim(40*[-1, 1])
    % ylim(5*[-1, 1]);
    % colorbar

    % Captura o frame e converte para imagem indexada
    frame1 = getframe(fig1);
    im1 = frame2im(frame1);
    [imind1, cm1] = rgb2ind(im1, 256);

    % Captura o frame e converte para imagem indexada
    frame2 = getframe(fig2);
    im2 = frame2im(frame2);
    [imind2, cm2] = rgb2ind(im2, 256);
    
    % Escreve no GIF
    if t == 1
        imwrite(imind1, cm1, filename1, 'gif', 'Loopcount', inf, 'DelayTime', delay_time);
        imwrite(imind2, cm2, filename2, 'gif', 'Loopcount', inf, 'DelayTime', delay_time);
    else
        imwrite(imind1, cm1, filename1, 'gif', 'WriteMode', 'append', 'DelayTime', delay_time);
        imwrite(imind2, cm2, filename2, 'gif', 'WriteMode', 'append', 'DelayTime', delay_time);
    end
end



%% Plot
execution_time = cputime - t0




