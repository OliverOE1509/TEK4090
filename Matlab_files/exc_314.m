clear all; close all; clc;

A = [9.9985 * 10^-1, 9.8510 * 10^-3;
    -2.9553 * 10^-2, 9.7030 * 10^-1];
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
function [x_hat_plus, P_plus] = runKalmanFilter(A, B, H, Q, R, y_meas, x_hat0_minus, P0_minus, N)
    x_hat_minus = zeros(2, N);
    x_hat_plus = zeros(2, N);
    P_minus = zeros(2,2,N);
    P_plus = zeros(2,2,N);

    x_hat_minus(:,1) = x_hat0_minus;
    P_minus(:, :, 1) = P0_minus;
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

[x_hat_plus_1, P_plus_1] = runKalmanFilter(A, B, H, Q, R, y_meas, [1.5; 0.5], [0.0001 0; 0 0.0001], N);
[x_hat_plus_2, P_plus_2] = runKalmanFilter(A, B, H, Q, R, y_meas, [50; -50], [100 0; 0 100], N);
[x_hat_plus_3, P_plus_3] = runKalmanFilter(A, B, H, Q, R, y_meas, [2; -1], [0.5 0.1; 0.1 0.5], N);
[x_hat_plus_4, P_plus_4] = runKalmanFilter(A, B, H, Q, R, y_meas, [1; 1], [0.001 0; 0 0.001], N);

error_1 = x_true - x_hat_plus_1;
error_2 = x_true - x_hat_plus_2;
error_3 = x_true - x_hat_plus_3;
error_4 = x_true - x_hat_plus_4;

dist1 = sqrt(sum(error_1.^2, 1));
dist2 = sqrt(sum(error_2.^2, 1));
dist3 = sqrt(sum(error_3.^2, 1));
dist4 = sqrt(sum(error_4.^2, 1));

figure;
plot(time, dist1, 'r-', 'LineWidth', 1.5); hold on;
plot(time, dist2, 'g-', 'LineWidth', 1.5);
plot(time, dist3, 'b-', 'LineWidth', 1.5);
plot(time, dist4, 'm-', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Euclidean Distance');
title('Distance Between True and Estimated State - Interesting Scenarios');
legend(...
    'Case 1: Slightly Wrong + Overconfident', ...
    'Case 2: Very Wrong + Very Uncertain', ...
    'Case 3: Somewhat Wrong + Realistic', ...
    'Case 4: Perfect + Overconfident', ...
    'Location', 'best');
grid on;

figure('Position', [80, 80, 900, 1200]);

subplot(4,2,1);
sigma_bounds_1_state1 = 3 * squeeze(sqrt(P_plus_1(1,1,:)));
plot(time, error_1(1,:), 'b', 'LineWidth', 1.5); hold on;
plot(time, sigma_bounds_1_state1, 'r--', 'LineWidth', 1);
plot(time, -sigma_bounds_1_state1, 'r--', 'LineWidth', 1);
ylabel('Error');
title('Case 1: State 1 Error (x0=[1.5; 0.5], P=[0.0001 0; 0 0.0001])');
legend('Error', '3σ Bound', 'Location', 'best');
grid on;

subplot(4,2,2);
sigma_bounds_1_state2 = 3 * squeeze(sqrt(P_plus_1(2,2,:)));
plot(time, error_1(2,:), 'b', 'LineWidth', 1.5); hold on;
plot(time, sigma_bounds_1_state2, 'r--', 'LineWidth', 1);
plot(time, -sigma_bounds_1_state2, 'r--', 'LineWidth', 1);
ylabel('Error');
title('Case 1: State 2 Error (x0=[1.5; 0.5], P=[0.0001 0; 0 0.0001])');
legend('Error', '3σ Bound', 'Location', 'best');
grid on;

subplot(4,2,3);
sigma_bounds_2_state1 = 3 * squeeze(sqrt(P_plus_2(1,1,:)));
plot(time, error_2(1,:), 'b', 'LineWidth', 1.5); hold on;
plot(time, sigma_bounds_2_state1, 'r--', 'LineWidth', 1);
plot(time, -sigma_bounds_2_state1, 'r--', 'LineWidth', 1);
ylabel('Error');
title('Case 2: State 1 Error (x0=[50; -50], P=[100 0; 0 100])');
legend('Error', '3σ Bound', 'Location', 'best');
grid on;

subplot(4,2,4);
sigma_bounds_2_state2 = 3 * squeeze(sqrt(P_plus_2(2,2,:)));
plot(time, error_2(2,:), 'b', 'LineWidth', 1.5); hold on;
plot(time, sigma_bounds_2_state2, 'r--', 'LineWidth', 1);
plot(time, -sigma_bounds_2_state2, 'r--', 'LineWidth', 1);
ylabel('Error');
title('Case 2: State 2 Error (x0=[50; -50], P=[100 0; 0 100])');
legend('Error', '3σ Bound', 'Location', 'best');
grid on;

subplot(4,2,5);
sigma_bounds_3_state1 = 3 * squeeze(sqrt(P_plus_3(1,1,:)));
plot(time, error_3(1,:), 'b', 'LineWidth', 1.5); hold on;
plot(time, sigma_bounds_3_state1, 'r--', 'LineWidth', 1);
plot(time, -sigma_bounds_3_state1, 'r--', 'LineWidth', 1);
ylabel('Error');
title('Case 3: State 1 Error (x0=[2; -1], P=[0.5 0.1; 0.1 0.5])');
legend('Error', '3σ Bound', 'Location', 'best');
grid on;


subplot(4,2,6);
sigma_bounds_3_state2 = 3 * squeeze(sqrt(P_plus_3(2,2,:)));
plot(time, error_3(2,:), 'b', 'LineWidth', 1.5); hold on;
plot(time, sigma_bounds_3_state2, 'r--', 'LineWidth', 1);
plot(time, -sigma_bounds_3_state2, 'r--', 'LineWidth', 1);
ylabel('Error');
title('Case 3: State 2 Error (x0=[2; -1], P=[0.5 0.1; 0.1 0.5])');
legend('Error', '3σ Bound', 'Location', 'best');
grid on;

subplot(4,2,7);
sigma_bounds_4_state1 = 3 * squeeze(sqrt(P_plus_4(1,1,:)));
plot(time, error_4(1,:), 'b', 'LineWidth', 1.5); hold on;
plot(time, sigma_bounds_4_state1, 'r--', 'LineWidth', 1);
plot(time, -sigma_bounds_4_state1, 'r--', 'LineWidth', 1);
ylabel('Error');
title('Case 4: State 1 Error (x0=[1; 1], P=[0.001 0; 0 0.001])');
legend('Error', '3σ Bound', 'Location', 'best');
grid on;

subplot(4,2,8);
sigma_bounds_4_state2 = 3 * squeeze(sqrt(P_plus_4(2,2,:)));
plot(time, error_4(2,:), 'b', 'LineWidth', 1.5); hold on;
plot(time, sigma_bounds_4_state2, 'r--', 'LineWidth', 1);
plot(time, -sigma_bounds_4_state2, 'r--', 'LineWidth', 1);
ylabel('Error');
title('Case 4: State 2 Error (x0=[1; 1], P=[0.001 0; 0 0.001])');
legend('Error', '3σ Bound', 'Location', 'best');
grid on;
