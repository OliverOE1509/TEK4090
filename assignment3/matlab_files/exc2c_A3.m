sigma_good_voltmeter = 1;
sigma_bad_voltmeter = 16;
N = 100;
R0 =  1250; % Ohms

function R_hat_WLS = generate_Rhat_WLS(N, R0, sigma_good, sigma_bad)
    i0 = -0.01 + 0.02 * rand(N, 1);    
    nu1 = sigma_good * randn(N, 1);
    nu2 = sigma_bad * randn(N, 1);
 
    % Correct: voltage = R0 * current + noise
    u_good = R0 * i0 + nu1; % Good voltmeter measurements
    u_bad = R0 * i0 + nu2;  % Bad voltmeter measurements
    
    % Use the KNOWN noise variances as weights (not computed from data)
    w_good = sigma_good^2;  % = 1
    w_bad = sigma_bad^2;    % = 256
    
    % Combine all measurements
    u_cat = [u_good; u_bad];
    i_cat = [i0; i0];  % Same current measurements repeated for both voltmeters
    
    numerator = 0;
    denominator = 0;
    
    for t = 1:(2*N)
        if t <= N
            weight = w_good;
        else
            weight = w_bad;
        end
        
        numerator = numerator + (i_cat(t) * u_cat(t)) / weight;
        denominator = denominator + (i_cat(t)^2) / weight;
    end
    
    R_hat_WLS = numerator / denominator;
end

function R_hat_LS = generate_Rhat_LS(N, R0, sigma_good, sigma_bad)
    i0 = -0.01 + 0.02 * rand(N, 1);    
    nu1 = sigma_good * randn(N, 1);
    nu2 = sigma_bad * randn(N, 1);
 
    % Correct: voltage = R0 * current + noise
    u_good = R0 * i0 + nu1; % Good voltmeter measurements
    u_bad = R0 * i0 + nu2;  % Bad voltmeter measurements
    
    % Combine all measurements
    u_cat = [u_good; u_bad];
    i_cat = [i0; i0];  % Same current measurements repeated for both voltmeters
    
    numerator = 0;
    denominator = 0;
    
    for t = 1:(2*N)
        numerator = numerator + (i_cat(t) * u_cat(t));
        denominator = denominator + (i_cat(t)^2);
    end
    
    R_hat_LS = numerator / denominator;
end


function [var_WLS_theoretical, var_LS_theoretical] = calculate_theoretical_variances(N)
    % Generate one set of current measurements
    i0 = -0.01 + 0.02 * rand(N, 1);
    i_cat = [i0; i0];  % 2N measurements
    
    % For WLS: Var(Ř_WLS) = 1 / sum_{t=1}^{2N} [i(t)^2 / w(t)]
    denominator_WLS = 0;
    for t = 1:(2*N)
        if t <= N
            weight = 1;  % w(t) = 1 for good voltmeter
        else
            weight = 256; % w(t) = 256 for bad voltmeter
        end
        denominator_WLS = denominator_WLS + (i_cat(t)^2) / weight;
    end
    var_WLS_theoretical = 1 / denominator_WLS;
    
    % For LS: Var(Ř_LS) = sum_{t=1}^{2N} [i(t)^2 σ_u(t)^2] / (sum_{t=1}^{2N} i(t)^2)^2
    numerator_LS = 0;
    denominator_LS = 0;
    for t = 1:(2*N)
        if t <= N
            sigma_u2 = 1;  % σ_u²(t) = 1 for good voltmeter
        else
            sigma_u2 = 16; % σ_u²(t) = 16 for bad voltmeter
        end
        numerator_LS = numerator_LS + (i_cat(t)^2) * sigma_u2;
        denominator_LS = denominator_LS + i_cat(t)^2;
    end
    denominator_LS = denominator_LS^2;
    var_LS_theoretical = numerator_LS / denominator_LS;
end


R_WLS = zeros(1e5, 1);
R_LS = zeros(1e5, 1);
for i = 1:1e6
    R_WLS(i) = generate_Rhat_WLS(N, R0, sigma_good_voltmeter, sigma_bad_voltmeter);
    R_LS(i) = generate_Rhat_LS(N, R0, sigma_good_voltmeter, sigma_bad_voltmeter);
end


disp(var(R_WLS))
disp(var(R_LS))
[var_R_WLS, var_R_LS] = calculate_theoretical_variances(N);
fprintf('WLS variance')
disp(var_R_WLS)

fprintf('LS variance')
disp(var_R_LS)



return
figure;
h1 = histogram(R_WLS, 'FaceAlpha', 0.5, 'DisplayName', 'WLS');
hold on;
h2 = histogram(R_LS, 'FaceAlpha', 0.5, 'DisplayName', 'LS');
xlabel('Resistance (Ohms)');
ylabel('Frequency');
title('Comparison of R_{LS} and R_{WLS} Estimates');
legend;
grid on;