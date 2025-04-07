function [rp, tp, ap] = fresnel_coefficients_p(theta_i, n1, n2, n3, dk, lambda0)

qk = @(n1, nk) sqrt(nk^2 - n1^2*sin(theta_i).^2)/nk^2;

q1 = qk(n1, n1);

q2 = qk(n1, n2);
q3 = qk(n1, n3);

beta_2 = (2*pi*dk/lambda0)*sqrt(n2^2 - n1^2*sin(theta_i).^2);

m11 = cos(beta_2);
m12 = -1j*sin(beta_2)./q2;
m21 = -1j*q2.*sin(beta_2);
m22 = cos(beta_2);

a = (m11 + m12.*q3);
b = (m21 + m22.*q3);

c = a.*q1 + b;
d = a.*q1 - b;

rp = d./c;
tp = 2*q1./c;
ap = 4*q1.*real(a.*conj(b) - q3)./(c.*conj(c));

end
