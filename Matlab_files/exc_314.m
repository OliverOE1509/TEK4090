clear all; close all; clc;

A = [9.9985 * 10^-1, 9.8510 * 10^-3;
    -2.9553 * 10-2, 9.7030 * 10^-1];
B = [4.9502 * 10^-5;
    9.8510 * 10^-3];

H = [1, 0];

Q = 1;
R = 0.01;
N = 1001;
dt = 0.01;
time = (0:N-1)' * dt;

x0_true = [1; 1];

rng(42)
w = sqrt(Q) * randn(N,1);
v = sqrt(R) * randn(N,1);
x_true = zeros(2, N);
x_true(:,1) = x0_true;

for k = 1:N-1
    x_true(:,k+1) = A * x_true(:,k) + B * w(k);
end

y_meas = H * x_true + v';
function [x_hat_plus, P_plus, x_hat_minus, P_minus] = ...
    runKalmanFilter(A, B, H, Q, R, y_meas, x_hat0_minus, P0_minus, N)
    x_hat_minus = zeros(2, N);
    x_hat_plus = zeros(2, N);
    P_minus = zeros(2,2,N);
    P_plus = zeros(2,2,N);

    x_hat_minus(:,1) = x_hat0_minus;
    P0_minus(:, :, 1) = P0_minus;
    for k=1:N
        Kk = P_minus(:, :, k) * H' / (H * P_minus(:, :, k) * H' + R);
        x_hat_plus(:, k) = x_hat_minus(:, k) + Kk * (y_meas(k) - H * x_hat_minus(:, k));
        P_plus(:, :, k) = (eye(2) - Kk * H) * P_minus(:, :, k);

        if k < N
            x_hat_minus(:, k+1) = A * x_hat_plus(:, k);
            P_minus(:, :, k+1) = A * P_plus(:, :, k) * A' + B * Q * B';
        end
    end
end


% Case 1: Initial state [0,0], P0 = I
fprintf('Running Case 1...\n');
x_hat0_1 = [0; 0];
P0_1 = [1 1;
    1 1];
[x_hat_plus_1, P_plus_1] = runKalmanFilter(A, B, H, Q, R, y_meas, [0;0], [0 0;0 0], N);
[x_hat_plus_2, P_plus_2] = runKalmanFilter(A, B, H, Q, R, y_meas, [1;1], [0.01 0.01;0.01 0.01], N);
[x_hat_plus_3, P_plus_3] = runKalmanFilter(A, B, H, Q, R, y_meas, [100; 100], [1 1; 1 1], N);
error_1 = x_true - x_hat_plus_1;
error_2 = x_true - x_hat_plus_2;
error_3 = x_true - x_hat_plus_3;
dist1 = sqrt(sum((error_1).^2));
dist2 = sqrt(sum((error_2).^2));
dist3 = sqrt(sum((error_3).^2));
figure;
plot(1:N, [dist1.', dist2.', dist3.']); legend('est1', 'est2', 'est3')
xlim([1, N]);

return
sigma_bounds_1 = 3 * squeeze(sqrt([P_plus_1(1,1,:), P_plus_1(2,2,:)]))';

% Plot Case 1
figure('Position', [100, 100, 1200, 800]);

subplot(2,2,1);
plot(time, error_1(1,:), 'b', 'LineWidth', 1.5); hold on;
plot(time, sigma_bounds_1(:,1), 'r--', 'LineWidth', 1);
plot(time, -sigma_bounds_1(:,1), 'r--', 'LineWidth', 1);
xlabel('Time (s)'); ylabel('Error');
title('Case 1: State 1 Error (x\_hat0 = [0;0], P0 = I)');
legend('Error', '3σ Bound', 'Location', 'best');
grid on;

subplot(2,2,2);
plot(time, error_1(2,:), 'b', 'LineWidth', 1.5); hold on;
plot(time, sigma_bounds_1(:,2), 'r--', 'LineWidth', 1);
plot(time, -sigma_bounds_1(:,2), 'r--', 'LineWidth', 1);
xlabel('Time (s)'); ylabel('Error');
title('Case 1: State 2 Error (x\_hat0 = [0;0], P0 = I)');
legend('Error', '3σ Bound', 'Location', 'best');
grid on;


sgtitle('Kalman Filter Estimation Errors with 3σ Bounds');



sgtitle('Kalman Filter Estimation Errors with 3σ Bounds (Continued)');

%% Convergence Analysis
fprintf('\n=== Convergence Analysis ===\n');

% Check if errors remain within 3σ bounds
for i = 1:4
    switch i
        case 1
            error = error_1;
            sigma_bounds = sigma_bounds_1;
            case_name = 'Case 1';
    end
    
    % Check if errors are within 3σ bounds (after initial transient)
    steady_state_start = 501; % Start checking from middle of simulation
    within_bounds_state1 = sum(abs(error(1, steady_state_start:end)) <= sigma_bounds(steady_state_start:end, 1)');
    
    percentage_state1 = within_bounds_state1 / (N - steady_state_start + 1) * 100;
    
    fprintf('%s:\n', case_name);
    fprintf('  State 1: %.2f%% of errors within 3σ bounds\n', percentage_state1);
    
    % Check final error magnitude
    final_error_state1 = abs(error(1, end));
    fprintf('  Final errors: State1 = %.4f', final_error_state1);
end

%% Display System Information
fprintf('=== System Information ===\n');
fprintf('Sampling interval: %.3f s\n', dt);
fprintf('Simulation duration: %.1f s\n', time(end));
fprintf('Process noise variance Q: %.2f\n', Q);
fprintf('Measurement noise variance R: %.4f\n', R);
fprintf('Eigenvalues of Phi: %.4f, %.4f\n', eig(A));

