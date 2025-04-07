function E = tmm_3p_TM(E_TM, fx, z, n, d, lambda0)
% Transfer Matrix Method for a 3phase system, TM mode

lambda = lambda0./n;
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

H_as = @(z, lambda) exp(1j*2*pi*z*sqrt(1/lambda^2 - fx.^2));

E1 = -q(1,:).*(E_TM.*H_as(z1, lambda(1)) - rp.*E_TM.*H_as(-z1, lambda(1)));
E2 = -q(1,:).*(E_TM.*H_as(z2, lambda(2)) + rp.*E_TM.*H_as(-z2, lambda(2)));
E3 = -q(end,:).*tp.*E_TM.*H_as(z3, lambda(3));

E = [E1;E2;E3];

end