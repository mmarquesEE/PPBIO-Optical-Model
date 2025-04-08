function E = tmm_3p_TM(E_TM, fx, z, n, d, lambda0)
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