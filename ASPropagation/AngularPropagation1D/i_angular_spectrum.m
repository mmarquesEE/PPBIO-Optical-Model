function [U, x] = i_angular_spectrum(A, fx)

dx = 1/mean(diff(fx));
L = 2^nextpow2(length(fx));

U = ifft(A*length(fx), L, 2);
x = (dx/L)*(-L/2:L/2-1);

end
