function delta_neff = calculate_sensorgram_from_formula(n_analyte, n_bulk, n_metal_complex, d_analyte, lambda)
% Calculates the sensorgram response based on the specific analytical
% formula provided by the user in equations (6) and (7).

    % --- Define terms based on your formula's notation ---
    % Epsilon_2r: Real part of the metal's dielectric constant (layer 2 in a 3-layer model)
    epsilon_2r = real(n_metal_complex^2);
    
    % N3: Refractive index of the analyte layer
    N3_sq = n_analyte^2;
    
    % N4: Refractive index of the bulk/environmental medium
    N4 = n_bulk;
    N4_sq = N4^2;

    % d3: Thickness of the layer causing the change.
    % Based on the physics, we will interpret this as the analyte layer thickness, d2.
    d3 = d_analyte;

    % --- Calculate Equation (6) ---
    term1 = (2 * pi * d3) / lambda;
    
    % Note: (-epsilon_2r * N4^2) will be positive since epsilon_2r for gold is negative.
    term2_numerator = (-epsilon_2r * N4_sq)^(3/2);
    term2_denominator = (epsilon_2r - N4_sq)^2;
    term2 = term2_numerator / term2_denominator;

    term3_numerator = N3_sq - N4_sq;
    term3_denominator = N3_sq; % As per your formula
    term3 = term3_numerator / term3_denominator;

    N_s_eff = term1 * term2 * term3 + N4_sq;
    
    % --- Calculate Equation (7) ---
    delta_neff = N_s_eff - N4;
end
