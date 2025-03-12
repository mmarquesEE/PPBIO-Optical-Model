function [k,P] = gen_rays(R,O,rx,ry,n,type)

if string(type) == "full"
    theta = 2 * pi * [0, rand(1, n - 1)];
    r = [0, sqrt(rand(1, n - 1))];
else
    theta = linspace(0, 2*pi, n);
    r = 1;
end

x = rx * r .* cos(theta);
y = ry * r .* sin(theta);
z = zeros(1, n);

k_ = [zeros(1,n);zeros(1,n);ones(1,n)];
P_ = [x; y; z];

k  = R*k_;
P  = R*P_ + O;

end

