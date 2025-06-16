% --- NEW HELPER FUNCTION to Calculate Validation Error Metrics ---
function [offset, avg_abs_error] = calculate_validation_error(y_true, y_hat)
% Calculates the offset and offset-corrected average absolute error.
% y_true: The original, ground-truth data vector.
% y_hat: The estimated/reconstructed data vector.

    % Ensure inputs are column vectors for consistent calculation
    y_true = y_true(:);
    y_hat = y_hat(:);

    % Equation 1: Calculate the mean offset 'c'
    offset = mean(y_hat - y_true);
    
    % Equation 2: Calculate the offset-corrected average absolute error 'E_A'
    avg_abs_error = mean(abs(y_true - (y_hat - offset)));
end