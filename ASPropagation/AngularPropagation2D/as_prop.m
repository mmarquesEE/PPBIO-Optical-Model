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

theta0 = 5/180*pi;

k = 2*pi/lambda1;

w0 = lambda1/(pi*theta0);
z0 = lambda1/(pi*theta0^2);

%% Angular Spectrum Propagation

d_alpha = 0.02/180*pi;
ds = lambda0/d_alpha;

dz = d_tf/20;

s = ds*linspace(-.5,.5,2^14);

%%
theta = 68.8/180*pi;

x = s*cos(theta);
z = -s*sin(theta);

U = gb_complex_amplitude(0.3, x, z, k, w0, z0);

[Ai1, fs] = angular_spectrum(U, s);

%% Using R-TMM
zt = ((-40e-6):dz:10e-6)';
At = tmm_3p_TM(Ai1, fs, zt, [n1;n2;n3], d_tf, lambda0);

[Ut, sp] = i_angular_spectrum(At, fs);

%% Plot
execution_time = cputime - t0

figure(Position=[150,150,1500,300])
hold on
imagesc(-1e6*sp, 1e6*zt, abs(Ut))
plot(1e6*[min(sp),max(sp)], [0,0], 'r', 1e6*[min(sp),max(sp)], 1e6*d_tf*[1,1], 'r');
hold off

title(['$|E_{xz}|$ ($\lambda=',num2str(lambda0_um,2),...
    '\mu$m, $\theta = ',num2str(theta/pi*180,2),'$deg, ',...
    '$\theta_{||}=',num2str(theta0/pi*180),'$deg)'],...
    Interpreter='latex', FontSize=16)
xlabel('x($\mu$m)', Interpreter='latex', FontSize=16)
ylabel('z($\mu$m)', Interpreter='latex', FontSize=16);
axis equal
xlim(100*[-1, 1])
% ylim(5*[-1, 1]);
% colorbar


