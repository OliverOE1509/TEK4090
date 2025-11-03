function R_hat = wls_estimate(i, u, w)
    % WLS_ESTIMATE - Weighted Least Squares resistance estimation
    %
    % Inputs:
    %   i - current measurements (vector)
    %   u - voltage measurements (vector)  
    %   w - weights (variances, vector)
    %
    % Output:
    %   R_hat - estimated resistance
    disp(size(w))
    
    % Calculate numerator: sum(i*u / w)
    numerator = sum((i .* u) ./ w);
    
    % Calculate denominator: sum(i^2 / w)
    denominator = sum((i.^2) ./ w);
    
    % WLS estimate
    R_hat = numerator / denominator;
end

function [R_wls, R_ls, estimates_wls, estimates_ls] = resistance_estimation_experiment(N, num_experiments)

    
    R0 = 1250;  % True resistance
    
    % Preallocate arrays
    estimates_wls = zeros(num_experiments, 1);
    estimates_ls = zeros(num_experiments, 1);
    
    for exp_idx = 1:num_experiments
        % Generate true current (uniform distribution)
        i0 = 0.02 * rand(2*N, 1) - 0.01;  % Uniform [-0.01, 0.01] A
        
        % Generate true voltage (Ohm's law)
        u0 = R0 * i0;
        
        % Add measurement noise
        % First N samples: good voltmeter (σ² = 1)
        u_measured = u0;
        u_measured(1:N) = u_measured(1:N) + randn(N, 1);
        
        % Next N samples: bad voltmeter (σ² = 16)  
        u_measured(N+1:end) = u_measured(N+1:end) + 4 * randn(N, 1);
        
        % Define weights (variances)
        w = ones(2*N, 1);
        w(1:N) = 1;      % Good voltmeter weights
        w(N+1:end) = 16; % Bad voltmeter weights
        
        % WLS estimation
        estimates_wls(exp_idx) = wls_estimate(i0, u_measured, w);
        
        % LS estimation (all weights = 1)
        estimates_ls(exp_idx) = wls_estimate(i0, u_measured, ones(2*N, 1));
    end
    
    % Calculate mean estimates
    R_wls = mean(estimates_wls);
    R_ls = mean(estimates_ls);
end

% Run the experiment
N = 100;
num_experiments = 1e5;
[R_wls, R_ls, estimates_wls, estimates_ls] = resistance_estimation_experiment(N, num_experiments);



% Plot histograms

disp(size(estimates_ls))
disp(size(estimates_wls))

figure;


normplot(estimates_wls)
title('Normal probability plot for the WLS estimates')

return
h1 = histogram(estimates_wls, 'FaceAlpha', 0.5, 'DisplayName', 'WLS', BinMethod='scott');
hold on;
h2 = histogram(estimates_ls, 'FaceAlpha', 0.5, 'DisplayName', 'LS', BinMethod='scott');
xlabel('Resistance (Ohms)');
ylabel('Frequency');
title('Comparison of R_{LS} and R_{WLS} Estimates');
legend;
grid on;
return
