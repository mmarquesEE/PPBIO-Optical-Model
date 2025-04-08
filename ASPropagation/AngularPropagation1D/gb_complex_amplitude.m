function U = gb_complex_amplitude(P, x, z, k, w0, z0)

w = w0*sqrt(1 + (z/z0).^2);

if z ~= 0
    R = z.*(1 + (z0./z).^2);
    eta = atan(z/z0);
    
    phi = k*z + k*x.^2./(2*R) - eta;
else
    phi = 0;
end

U = sqrt(2*P/pi)./w.*exp(-x.^2./w.^2).*exp(-1j*phi);

end