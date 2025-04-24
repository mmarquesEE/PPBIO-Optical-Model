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
% n3 = 1;
dn3 = 0.02*sin(2*pi*linspace(0,1,100))-0.015;

% thickness of gold layer
d_tf = 50e-9;

k = 2*pi/lambda0;

%% Utils
proj = @(ki, theta) cos(theta)*ki + sin(theta)*sqrt(k^2 - ki.^2);
J = @(ki, theta) cos(theta) - sin(theta)*ki./sqrt(k^2 - ki.^2);
l_tilt = @(i,theta) -sin(theta)*i;

%% PPBIO
h = 3e-3; alpha = 34/180*pi;
d1 = 5e-3;
d2 = h/2*(1 + 3/cos(2*alpha));
L = 3e-3;

f = d1 + real(n1)*d2;

theta = 2*alpha;

%% Angular Spectrum - y direction
Ny = 4;
ky = k/100;
y_lim = 1.2*L/2;
y = linspace(-y_lim,y_lim,Ny);

arg_tay = y/L;
tay = abs(arg_tay) < 1/2;
tay(abs(arg_tay) == 1/2) = 1/2;

Ary = (tay*exp(-1j*transpose(y)*ky)).*exp(1j*k*2*f).*exp(-1j*pi*lambda0*2*f*(ky/2/pi).^2);
Ery = (Ary*exp(1j*transpose(ky)*y));
Ery = Ery/max(abs(Ery));

%% Angular Spectrum -  x direction

Nx = 2048;

x_lim = 1.2*L/2;
x = linspace(-x_lim,x_lim,Nx);

kx_lim = 0.1*k;
kx = linspace(-kx_lim,kx_lim,Nx);
ku = proj(kx, theta);

arg_tax = f*kx/k/L;
tax = abs(arg_tax) < 1/2;
tax(arg_tax == 1/2) = 1/2;

xi1 = sqrt(n1^2 - n1^2*(ku/k).^2);
xi2 = sqrt(n2^2 - n1^2*(ku/k).^2);

q1 = xi1/n1^2;
q2 = xi2/n2^2;

beta2 = k*d_tf*xi2;

m11 = cos(beta2);
m22 = m11;
m12 = -1j./q2.*sin(beta2);
m21 = -1j*q2.*sin(beta2);

Arp = sqrt(lambda0*f)*tax.*J(ku, theta).*J(kx,-theta)...
    .*exp(1j*k*f).*exp(-1j*pi*lambda0*f*(kx/2/pi).^2);

Erx2 = (Arp*exp( 1j*transpose(kx)*x));

[Ey2, Ex2] = meshgrid(Ery, Erx2);
I2 = abs(Ex2).*abs(Ey2);

filename = "test.gif"; delay_time = 0.1;

figure(Position=[100,100,1200,600]);
tiledlayout(2,2);

for i=1:length(dn3)
    xi3 = sqrt((n3+dn3(i))^2 - n1^2*(ku/k).^2);
    q3 = xi3/(n3+dn3(i))^2;

    rp = ((m11 + m12.*q3).*q1 - (m21 + m22.*q3))./...
         ((m11 + m12.*q3).*q1 + (m21 + m22.*q3));

    Arx1 = Arp.*rp;

    Erx1 = (Arx1*exp( 1j*transpose(kx)*x));

    [Ey1, Ex1] = meshgrid(Ery, Erx1);
    
    I1 = abs(Ex1).*abs(Ey1);

    Ir = I1./I2;
    Ir = Ir/max(Ir);
    
    fontsize = 16;
    
    nexttile(1,[2, 1]);
    imagesc(1e3*y, -1e3*x, Ir);
    colormap("gray");
    xlabel('y(mm)', Interpreter='latex', FontSize=fontsize);
    ylabel('x(mm)', Interpreter='latex', FontSize=fontsize);
    title('$\frac{|E_{analyte}|^2}{|E_{air}|^2}$', Interpreter='latex',FontSize=fontsize);
    axis equal
    xlim(1.5*[-1, 1]);
    ylim(1.5*[-1, 1]);
    
    nexttile(2);
    plot( ...
        -1e3*x, (abs(Erx1)/max(abs(Erx1))).^2, ...
        -1e3*x, (abs(Erx2)/max(abs(Erx2))).^2 ...
    );
    legend(["analyte","air"],Interpreter="latex");
    grid on
    xlabel('x(mm)', Interpreter='latex', FontSize=fontsize);
    ylabel('$|E|^2$',Interpreter='latex', FontSize=fontsize);
    title('$|E_{analyte}|^2$ for $(y=0)$ (normalized)', Interpreter='latex',FontSize=fontsize);

    nexttile(4);
    plot(1e3*dn3(1:i));
    xlim([1, length(dn3)]);
    ylim(1e3*[min(dn3),max(dn3)])
    grid on;
    xlabel('samples',Interpreter='latex',FontSize=fontsize);
    ylabel('$10^3\times\Delta n_{analyte}$',Interpreter='latex',FontSize=fontsize);
    title( ...
        "Analyte's refractive index change ($n_{water}="+num2str(real(n3),3)+"$)", ...
        Interpreter="latex", ...
        FontSize=fontsize ...
    );
    

    % Captura o frame e converte para imagem indexada
    frame = getframe(gcf);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    
    % Escreve no GIF
    if i == 1
        imwrite(imind, cm, filename, 'gif', 'Loopcount', inf, 'DelayTime', delay_time);
    else
        imwrite(imind, cm, filename, 'gif', 'WriteMode', 'append', 'DelayTime', delay_time);
    end
end

%%
elapsed_time = cputime - t0
