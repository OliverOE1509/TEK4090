% ========== PART 1(i): PARAMETER UNCERTAINTY - FEEDFORWARD ONLY ==========

% Original nominal mass
M_nominal = 10;

% Define uncertainty cases
uncertainty_cases = {
    'Nominal', M_nominal, M_nominal;                    % No uncertainty
    'M1 +1%',  M_nominal*1.01, M_nominal;              % M1 increased 1%
    'M1 -1%',  M_nominal*0.99, M_nominal;              % M1 decreased 1%  
    'M2 +1%',  M_nominal, M_nominal*1.01;              % M2 increased 1%
    'M2 -1%',  M_nominal, M_nominal*0.99;              % M2 decreased 1%
    'Both +1%', M_nominal*1.01, M_nominal*1.01;        % Both increased 1%
    'Both -1%', M_nominal*0.99, M_nominal*0.99         % Both decreased 1%
};

% Store results for comparison
results = struct();



for case_idx = 1:size(uncertainty_cases, 1)
    case_name = uncertainty_cases{case_idx, 1};
    M1_actual = uncertainty_cases{case_idx, 2};
    M2_actual = uncertainty_cases{case_idx, 3};
    
    fprintf('Simulating case: %s (M1=%.2f, M2=%.2f)\n', case_name, M1_actual, M2_actual);
    
    % ========== SYSTEM WITH UNCERTAINTY ==========
    % PLANT (actual system with uncertainty)
    denom_actual = M1_actual*M2_actual + 2*(M1_actual + M2_actual) + 3;
    k_actual = (M1_actual + M2_actual + 3) / denom_actual;
    c_actual = 0.1 * (M1_actual + M2_actual + 3) / denom_actual;
    
    A_actual = [0, 0, 1, 0;
               0, 0, 0, 1;
               -k_actual, k_actual, -c_actual, c_actual;
               k_actual, -k_actual, c_actual, -c_actual];
    
    B_actual = [0; 0; (M2_actual + 2)/denom_actual; -1/denom_actual];
    
    % ========== SIMULATION WITH MODEL MISMATCH - FEEDFORWARD ONLY ==========
    x0 = [0; 0; 0; 0];
    x = x0;
    X_sim_actual = []; U_total_actual = []; Y_sim_actual = [];
    
    for k = 1:length(t)
        % FEEDFORWARD ONLY - no feedback component
        u_total = u_inv(k);  % ONLY feedforward, no feedback
        
        % Simulate with ACTUAL plant
        X_sim_actual = [X_sim_actual, x];
        U_total_actual = [U_total_actual, u_total];
        Y_sim_actual = [Y_sim_actual, x(2)];
        
        if k < length(t)
            dt = t(k+1) - t(k);
            xdot = A_actual * x + B_actual * u_total;
            x = x + xdot * dt;
        end
    end
    
    % Store results
    tracking_error_actual = y_d - Y_sim_actual;
    results(case_idx).name = case_name;
    results(case_idx).RMSE = sqrt(mean(tracking_error_actual.^2));
    results(case_idx).max_error = max(abs(tracking_error_actual));
    results(case_idx).control_effort = trapz(abs(U_total_actual));
    results(case_idx).Y_sim = Y_sim_actual;
    results(case_idx).U_total = U_total_actual;
    results(case_idx).tracking_error = tracking_error_actual;
end

% ========== 3-PLOT COMPARISON ==========
figure('Position', [100, 100, 1200, 600]);

% Plot 1: Output tracking trajectories
subplot(1, 3, 1);
plot(t, y_d, 'k--', 'LineWidth', 3, 'DisplayName', 'Desired'); hold on;
colors = lines(length(results));
for i = 1:length(results)
    plot(t, results(i).Y_sim, '-', 'LineWidth', 1.5, ...
         'Color', colors(i,:), 'DisplayName', results(i).name);
end
title('Output Tracking');
xlabel('Time (s)'); ylabel('Position');
legend('Location', 'best'); grid on;

% Plot 2: Control efforts
subplot(1, 3, 2);
for i = 1:length(results)
    plot(t, results(i).U_total, '-', 'LineWidth', 1.5, ...
         'Color', colors(i,:), 'DisplayName', results(i).name); hold on;
end
title('Control Efforts');
xlabel('Time (s)'); ylabel('Control Input');
legend('Location', 'best'); grid on;

% Plot 3: Tracking errors
subplot(1, 3, 3);
for i = 1:length(results)
    plot(t, results(i).tracking_error, '-', 'LineWidth', 1.5, ...
         'Color', colors(i,:), 'DisplayName', results(i).name); hold on;
end
title('Tracking Errors');
xlabel('Time (s)'); ylabel('Error');
legend('Location', 'best'); grid on;

sgtitle('Feedforward Only: Robustness to Parameter Uncertainty');
