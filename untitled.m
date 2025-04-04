clear all
close all
clc

% phase_grad = @(x, z) [k*x./R(z);k - k*x.^2./(2*R(z).^2).*(1 - z0^2./z.^2) - w0^2./w(z).^2];
% phase_grad_s = phase_grad(x, z);
% k_vec = phase_grad_s/vecnorm(phase_grad_s);

% ns = [sin(theta);cos(theta)];
% 
% theta_i = acos(transpose(ns)*k_vec);

%% 
lambda0 = 670e-9;

n = 1.51;
theta0 = 5/180*pi;
lambda = lambda0/n;

k = 2*pi/lambda;

w0 = lambda/(pi*theta0);
z0 = lambda/(pi*theta0^2);

%%

slim = 128*lambda;
s = slim*linspace(-1,1,4096);

theta = 68/180*pi;

x = s*cos(theta);
z = -s*sin(theta);

U = gb_complex_amplitude(x, z, k, w0, z0);

[A, fs] = angular_spectrum(U, s);
% ignore evanescent fields
alpha = fs*lambda;

A = A(abs(alpha) < 1);
ang = asin(alpha(abs(alpha) < 1));

[rp, tp, ~] = fresnel_coefficients_p();

figure
subplot(211);
plot(180*ang/pi,abs(A));
subplot(212);
plot(180*ang/pi,(180/pi)*angle(A));

%%
function U = gb_complex_amplitude(x, z, k, w0, z0)

w = w0*sqrt(1 + (z/z0).^2);

if z ~= 0
    R = z.*(1 + (z0./z).^2);
    eta = atan(z/z0);
    
    phi = k*z + k*x.^2./(2*R) - eta;
else
    phi = 0;
end

U = w0./w.*exp(-x.^2./w.^2).*exp(-1j*phi);

end

function [rp, tp, ap] = fresnel_coefficients_p(theta_i, n1, nk, dk, lambda0)

qk = @(n1, nk) sqrt(nk^2 - n1^2*sin(theta_i).^2)/nk^2;

q1 = qk(n1, n1);

qn = qk(n1, nk);
beta_k = (2*pi*dk/lambda0)*sqrt(nk^2 - n1^2*sin(theta_i).^2);

m11 = cos(beta_k);
m12 = -1j*sin(beta_k)/qn;
m21 = -1j*qn*sin(beta_k);
m22 = cos(beta_k);

a = (m11 + m12.*qn);
b = (m21 + m22.*qn);

c = a.*q1 + b;
d = a.*q1 - b;

rp = d./c;
tp = 2*q1./c;
ap = 4*q1.*real(a.*conj(b) - qn)./(c.*conj(c));

end

function [A, fx] = angular_spectrum(U, x)

Fs = 1/mean(diff(x));
L = length(U);

A  = fftshift(fft(U, L))/L;
fx = (Fs/L)*(-L/2:L/2-1);

end



