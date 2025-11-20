sigma_good_voltmeter = 1;
sigma_bad_voltmeter = 16;
N = 100;
R0 = 1250; % Ohms

function R_hat_WLS = generate_Rhat_WLS(N, R0, sigma_good, sigma_bad)
    i0 = -0.01 + 0.02 * rand(N, 1);    
    nu1 = sigma_good * randn(N, 1);
    nu2 = sigma_bad * randn(N, 1);
 
    % Correct: voltage = R0 * current + noise
    u_good = R0 * i0 + nu1; % Good voltmeter measurements
    u_bad = R0 * i0 + nu2;  % Bad voltmeter measurements
    
    % Use the KNOWN noise variances as weights
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
    
    % Check for denominator issues
    if denominator == 0
        R_hat_LS = R0; % Fallback to true value
    else
        R_hat_LS = numerator / denominator;
    end
end

function var_R_LS_theoretical = find_var_LS_theoretical(N)
    % Theoretical variance calculation using your derived formula:
    % Var(Ř_LS) = sum[i(t)^2 σ_u(t)^2] / (sum[i(t)^2])^2
    
    % For uniform distribution i ~ U[-0.01, 0.01]
    E_i2 = ((0.01 - (-0.01))^2) / 12;  % = 3.333e-5
    
    % Expected numerator: N*E[i^2]*1^2 + N*E[i^2]*16^2 = N*E_i2 + 256*N*E_i2
    expected_numerator = N * E_i2 * (1^2) + N * E_i2 * (16^2);
    expected_numerator = 257 * N * E_i2;
    
    % Expected denominator: (2*N*E_i2)^2
    expected_denominator = (2 * N * E_i2)^2;
    
    var_R_LS_theoretical = expected_numerator / expected_denominator;
end

function var_R_WLS_theoretical = find_var_WLS_theoretical(N)
    % Theoretical variance calculation using your derived formula:
    % Var(Ř_WLS) = 1 / sum[i(t)^2 / w(t)]
    
    % For uniform distribution i ~ U[-0.01, 0.01]
    E_i2 = ((0.01 - (-0.01))^2) / 12;  % = 3.333e-5
    
    % Expected denominator: N*E[i^2/1] + N*E[i^2/256] = N*E_i2 + N*E_i2/256
    expected_denominator = N * E_i2 + N * (E_i2 / 256);
    
    var_R_WLS_theoretical = 1 / expected_denominator;
end

% Run experiments with smaller sample size first to debug
num_experiments = 1000; % Start small for debugging
R_WLS = zeros(num_experiments, 1);
R_LS = zeros(num_experiments, 1);

fprintf('Running %d experiments...\n', num_experiments);
for i = 1:num_experiments
    R_WLS(i) = generate_Rhat_WLS(N, R0, sigma_good_voltmeter, sigma_bad_voltmeter);
    R_LS(i) = generate_Rhat_LS(N, R0, sigma_good_voltmeter, sigma_bad_voltmeter);
    
    % Debug: print first few values
    if i <= 5
        fprintf('Experiment %d: WLS = %.2f, LS = %.2f\n', i, R_WLS(i), R_LS(i));
    end
end

% Remove outliers (more than 3 standard deviations from mean)
LS_mean = mean(R_LS);
LS_std = std(R_LS);
valid_indices = abs(R_LS - LS_mean) < 3 * LS_std;

R_WLS_clean = R_WLS(valid_indices);
R_LS_clean = R_LS(valid_indices);

% Calculate experimental variances
var_WLS_experimental = var(R_WLS_clean);
var_LS_experimental = var(R_LS_clean);

% Calculate theoretical variances
var_WLS_theoretical = find_var_WLS_theoretical(N);
var_LS_theoretical = find_var_LS_theoretical(N);

% Display results
fprintf('\n=== RESULTS ===\n');
fprintf('Experimental Results (after removing outliers):\n');
fprintf('WLS variance: %.2f\n', var_WLS_experimental);
fprintf('LS variance: %.2f\n', var_LS_experimental);
fprintf('Original data points: %d, Clean data points: %d\n', num_experiments, sum(valid_indices));

fprintf('\nTheoretical Results:\n');
fprintf('WLS variance: %.2f\n', var_WLS_theoretical);
fprintf('LS variance: %.2f\n', var_LS_theoretical);

fprintf('\nComparison (Error %%):\n');
fprintf('WLS: %.1f%% error\n', 100*abs(var_WLS_theoretical - var_WLS_experimental)/var_WLS_theoretical);
fprintf('LS: %.1f%% error\n', 100*abs(var_LS_theoretical - var_LS_experimental)/var_LS_theoretical);

% Plot histograms
figure;
subplot(1,2,1);
histogram(R_WLS_clean, 50, 'FaceAlpha', 0.7, 'FaceColor', 'blue');
xlabel('Resistance (Ohms)');
ylabel('Frequency');
title(sprintf('WLS Estimates\nVar = %.1f', var_WLS_experimental));
grid on;
xlim([1200, 1300]);

subplot(1,2,2);
histogram(R_LS_clean, 50, 'FaceAlpha', 0.7, 'FaceColor', 'red');
xlabel('Resistance (Ohms)');
ylabel('Frequency');
title(sprintf('LS Estimates\nVar = %.1f', var_LS_experimental));
grid on;
xlim([1200, 1300]);

% Combined plot
figure;
histogram(R_WLS_clean, 50, 'FaceAlpha', 0.5, 'FaceColor', 'blue', 'DisplayName', sprintf('WLS (Var=%.1f)', var_WLS_experimental));
hold on;
histogram(R_LS_clean, 50, 'FaceAlpha', 0.5, 'FaceColor', 'red', 'DisplayName', sprintf('LS (Var=%.1f)', var_LS_experimental));
xlabel('Resistance (Ohms)');
ylabel('Frequency');
title('Comparison of R_{LS} and R_{WLS} Estimates');
legend;
grid on;
xlim([1200, 1300]);