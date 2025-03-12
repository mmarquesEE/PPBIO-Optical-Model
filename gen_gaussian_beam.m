function [k,P] = gen_gaussian_beam(R,O,t0x,t0y,lambda,na,np)

w0x = lambda/(pi*t0x);
w0y = lambda/(pi*t0y);

a_ = linspace(0.,2*pi,na);
p_ = linspace(0.,1.,np);

[A,B] = meshgrid(a_,p_);
a = reshape(A,1,[]); b = a + pi/4;
p = reshape(B,1,[]);

x0 = 0; y0 = 0; l0 = 0; m0 = 0;

x1 =  p .* w0x .* cos(a);
y1 =  p .* w0y .* sin(a);
l1 = -p .* t0x .* sin(a);
m1 =  p .* t0y .* cos(a);

x2 =  p .* w0x .* cos(a + pi/2);
y2 =  p .* w0y .* sin(a + pi/2);
l2 = -p .* t0x .* sin(a + pi/2);
m2 =  p .* t0y .* cos(a + pi/2);

x3 = -p .* w0x .* cos(b);
y3 = -p .* w0y .* sin(b);
l3 = -p .* t0x .* sin(b);
m3 =  p .* t0y .* cos(b);

x4 = -p .* w0x .* cos(b + pi/2);
y4 = -p .* w0y .* sin(b + pi/2);
l4 = -p .* t0x .* sin(b + pi/2);
m4 =  p .* t0y .* cos(b + pi/2);

position = @(x,y) [x;y;zeros(1,length(x))]; 

dir_cossines = @(l,m) [
    l./sqrt(l.^2 + m.^2 + 1);
    m./sqrt(l.^2 + m.^2 + 1);
    1./sqrt(l.^2 + m.^2 + 1);
];

P = R*position([x0,x1,x2,x3,x4],[y0,y1,y2,y3,y4]) + O;
k = R*dir_cossines([l0,l1,l2,l3,l4],[m0,m1,m2,m3,m4]);

end

