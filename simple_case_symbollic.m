clear all
clc

% Inclua o pacote de simbolos
syms lambda theta0 z x f_c f_f n_p theta s real

W0 = lambda/(pi*theta0);
z0 = lambda/(pi*theta0^2);

z0c = (f_c^2/(2*z0^2))*z0;
z0f = z0c/(1 + (z0c/f_f)^2);
z0r = n_p*z0f;

k = 2*pi*n_p/lambda;

W0r = sqrt(lambda*z0f/pi);

W = W0r*sqrt(1 + (z/z0r)^2);
R = z*(1 + (z0r/z)^2);
zetaa = atan(z/z0r);

varphi = simplify(k*z - zetaa + (k*x^2)/(2*R));
grad_varphi = simplify([(k*x)/R; k - (k*x^2)/(2*R^2)*(1 - (z0r/z)^2) - (W0r/W)^2])

latex(grad_varphi)
u = simplify((W0r/W)*exp(-(x/W)^2)*exp(-1j*varphi));

Varphi = simplify(subs(varphi, [x, z], [s*cos(theta), -s*sin(theta)]))
Grad_varphi = simplify(subs(grad_varphi, [x, z], [s*cos(theta), -s*sin(theta)]))
U = simplify(subs(u, [x, z], [s*cos(theta), -s*sin(theta)]))
