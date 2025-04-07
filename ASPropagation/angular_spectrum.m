function [A, fx] = angular_spectrum(U, x)

Fs = 1/mean(diff(x));
L = 2^nextpow2(length(x));

A  = fft(U, L, 2)/L;
fx = (Fs/L)*[0,(1:(L/2-1)),-(1:L/2)];

end