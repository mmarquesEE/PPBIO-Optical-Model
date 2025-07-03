function [Rp, resonance_angle] = fresnel_spr_curve(angles_deg, n0, n1, n2, n3, d1, d2, wavelength)
% Calculates the SPR reflectivity curve for p-polarized light in a 4-layer system.
% This is a MATLAB implementation of the provided Python/JAX Fresnel model.

    % Convert input angles from degrees to radians
    th = deg2rad(angles_deg);
    
    % Initialize output array for reflectivity
    Rp = zeros(size(th));
    
    % This loop is the equivalent of jax.vmap, applying the calculation for each angle
    parfor i = 1:length(th)
        current_th = th(i);
        
        % Check for TIR condition to avoid complex numbers in sqrt where not needed
        sin_th_sq = (n0 * sin(current_th))^2;
        
        % Compute q values for each layer (using complex numbers for generality)
        q0 = sqrt(n0^2 - sin_th_sq + 0i) / n0^2;
        q1 = sqrt(n1^2 - sin_th_sq + 0i) / n1^2;
        q2 = sqrt(n2^2 - sin_th_sq + 0i) / n2^2;
        q3 = sqrt(n3^2 - sin_th_sq + 0i) / n3^2;
        
        % Compute beta values for layers 1 and 2
        beta1 = 2 * pi * d1 * sqrt(n1^2 - sin_th_sq + 0i) / wavelength;
        beta2 = 2 * pi * d2 * sqrt(n2^2 - sin_th_sq + 0i) / wavelength;
        
        % Layer matrices M1 and M2
        M1 = [cos(beta1), -1j * sin(beta1) / q1; 
              -1j * q1 * sin(beta1), cos(beta1)];
          
        M2 = [cos(beta2), -1j * sin(beta2) / q2; 
              -1j * q2 * sin(beta2), cos(beta2)];
        
        % Overall matrix product
        M = M1 * M2;
        
        % Reflection coefficient for p-polarized light
        numerator = (M(1,1) + M(1,2) * q3) * q0 - (M(2,1) + M(2,2) * q3);
        denominator = (M(1,1) + M(1,2) * q3) * q0 + (M(2,1) + M(2,2) * q3);
        rp = numerator / denominator;
        
        % Reflectivity is the squared magnitude of the reflection coefficient
        Rp(i) = abs(rp)^2;
    end
    
    % Find the resonance angle (angle of minimum reflectivity)
    [~, min_idx] = min(Rp);
    resonance_angle = angles_deg(min_idx);
end