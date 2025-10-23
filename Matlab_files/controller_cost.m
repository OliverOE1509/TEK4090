%% LQR Cost Calculation - Similar to your screenshot
clear; clc; close all;

% System matrices (from your screenshot)
A = [-0.5 -1; 1 0];
B = [1; 0];

% LQR weights (from your screenshot)
Q = [1 0; 0 1];
R = 0.1;

% Calculate optimal gain using LQR
[K, P, E] = lqr(A, B, Q, R);

fprintf('Optimal gain K = [%.3f, %.3f]\n', K(1), K(2));
fprintf('Your K1 = %.3f\n', K(1));

% Calculate the optimal cost for a given initial condition
x0 = [10; 5];  % Example initial state

% Optimal cost J = x₀ᵀ P x₀
J_optimal = x0' * P * x0;

fprintf('Optimal cost J = x₀ᵀ P x₀ = %.1f\n', J_optimal);
fprintf('Your J = 637.3\n');

%% Display system information
fprintf('\n--- System Information ---\n');
disp('A = '); disp(A);
disp('B = '); disp(B);
disp('Q = '); disp(Q);
disp('R = '); disp(R);
disp('P = '); disp(P);

%% Simulate the closed-loop system
sys_open = ss(A, B, eye(2), 0);  % Original system
K_lqr = -K;  % Control law: u = -Kx

% Closed-loop system
A_cl = A + B * K_lqr;
sys_cl = ss(A_cl, B, eye(2), 0);

% Time vector
t = 0:0.1:20;

% Simulate response
[x, t] = initial(sys_cl, x0, t);

% Calculate control effort
u = zeros(length(t), 1);
for i = 1:length(t)
    u(i) = K_lqr * x(i,:)';
end

%% Plot results - Similar to your screenshot
figure('Position', [100, 100, 1200, 800]);

% Subplot 1: States
subplot(2,2,1);
plot(t, x(:,1), 'b-', 'LineWidth', 2); hold on;
plot(t, x(:,2), 'r-', 'LineWidth', 2);
xlabel('Time');
ylabel('State, x');
title('System States');
legend('x_1', 'x_2');
grid on;

% Subplot 2: Control Effort
subplot(2,2,2);
plot(t, u, 'g-', 'LineWidth', 2);
xlabel('Time');
ylabel('Controller Effort, u');
title('Control Effort');
grid on;

% Subplot 3: Phase portrait
subplot(2,2,3);
plot(x(:,1), x(:,2), 'k-', 'LineWidth', 2);
hold on;
plot(x0(1), x0(2), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'red');
plot(0, 0, 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'green');
xlabel('x_1');
ylabel('x_2');
title('Phase Portrait');
legend('Trajectory', 'Initial State', 'Equilibrium', 'Location', 'best');
grid on;

% Subplot 4: Cost components
subplot(2,2,4);
state_cost = zeros(length(t), 1);
control_cost = zeros(length(t), 1);
for i = 1:length(t)
    state_cost(i) = x(i,:) * Q * x(i,:)';
    control_cost(i) = u(i) * R * u(i)';
end
total_cost = state_cost + control_cost;

plot(t, state_cost, 'b-', 'LineWidth', 2); hold on;
plot(t, control_cost, 'r-', 'LineWidth', 2);
plot(t, total_cost, 'k--', 'LineWidth', 2);
xlabel('Time');
ylabel('Cost');
title('Cost Components');
legend('State Cost (x^TQx)', 'Control Cost (u^TRu)', 'Total Cost', 'Location', 'best');
grid on;

%% Calculate cumulative cost (verification)
% Discrete approximation of the integral
dt = t(2) - t(1);
J_cumulative = sum(total_cost) * dt;

fprintf('\n--- Cost Verification ---\n');
fprintf('Theoretical cost (x₀ᵀPx₀): %.4f\n', J_optimal);
fprintf('Simulated cumulative cost: %.4f\n', J_cumulative);
fprintf('Difference: %.4f\n', abs(J_optimal - J_cumulative));

%% Additional analysis: Eigenvalues
fprintf('\n--- Stability Analysis ---\n');
fprintf('Open-loop eigenvalues: %.3f ± %.3fi\n', real(eig(A)), imag(eig(A)));
fprintf('Closed-loop eigenvalues: %.3f ± %.3fi\n', real(eig(A_cl)), imag(eig(A_cl)));