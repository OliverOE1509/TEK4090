% Define parameters
num_P = [1, 1];
den_P = [2, -4, 13];
num_C = [1];
den_C = [1, 0];
K = 1;

% Call the function with separate arguments
isStable = check_stability_separate(num_P, den_P, num_C, den_C, K);

% Function definition for separate arguments
function isStable = check_stability_separate(num_P, den_P, num_C, den_C, K)
% CHECK_STABILITY_SEPARATE Check stability using separate arguments
    
    % Compute the convolutions
    den_product = conv(den_P, den_C);
    num_product = conv(num_P, num_C);
    
    % Ensure both polynomials have the same length by padding with zeros
    max_length = max(length(den_product), length(num_product));
    den_product_padded = [zeros(1, max_length - length(den_product)), den_product];
    num_product_padded = [zeros(1, max_length - length(num_product)), num_product];
    
    % Form the characteristic polynomial
    char_poly = den_product_padded + K * num_product_padded;
    
    % Find roots of the characteristic polynomial
    r = roots(char_poly);
    
    % Display the roots
    fprintf('For K = %.4f, roots are:\n', K);
    disp(r);
    
    % Check if all real parts are negative
    if all(real(r) < 0)
        fprintf('System is STABLE for K = %.4f.\n', K);
        isStable = true;
    else
        fprintf('System is UNSTABLE for K = %.4f.\n', K);
        isStable = false;
    end
end