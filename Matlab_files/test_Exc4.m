%% TEK4090 - Assignment 2, Exercise 4: Building Temperature Control

clear all; close all; clc;

%% Define model parameters
c = 1/(1.08*10^6);  % [K/J] - Inverse thermal capacitance

% Setpoint temperature 
x_bar = 20 + 273;  % [K] - Setpoint in Kelvin

% Initial conditions
x_ic = [16; 13; 18; 11; 20] + 273;  % [K]
z_ic = x_ic - x_bar;  % [K]

%% External input functions
d_t = @(t) max(0, 500*sin(2*pi*t/86400) - 300);
x_a_dev_t = @(t) 5*sin(2*pi*t/86400);

%% System matrices (FIXED: A matrix scaled by c)
A = [-30,  10,   0,   0,  10;
      10, -40,  10,   0,  10;
       0,  10, -40,  10,  10;
       0,   0,  10, -30,  10;
      10,  10,  10,  10, -50] * c;  % CRITICAL FIX!

B_u = [c, 0, 0;
       0, c, 0; 
       0, 0, 0;
       0, 0, c;
       0, 0, 0];

B_d = eye(5) * c;
B_a = [10; 10; 10; 10; 10] * c;
C = [1, 0, 0, 0, 0;
     0, 0, 1, 0, 0];
D = zeros(2, 3);

%% LQR Controller Design 
Q = eye(5) * 1e6;    % State penalty
R = eye(3) * 1;      % Control penalty

lqr_sys = ss(A, B_u, C, D);
[K, S, P] = lqr(lqr_sys, Q, R);

fprintf('LQR Design Complete:\n');
fprintf('Q = %.0e * I, R = %.0e * I\n', Q(1,1), R(1,1));

%% Observer Design
if rank(obsv(A, C)) == size(A,1)
    fprintf('System is observable.\n');
else
    error('System not observable - check C matrix');
end

% Use reasonable observer poles
desired_obs_poles = [-0.1, -0.2, -0.3, -0.4, -0.5];  % Slower poles
L = place(A', C', desired_obs_poles)';

%% Augmented System
A_aug = [A,          -B_u*K;
         L*C,   A - B_u*K - L*C];

B_ext_top = [B_d, B_a];
B_ext_bottom = [B_d, B_a];
B_aug = [B_ext_top; B_ext_bottom];

C_aug = [C, zeros(size(C))];
D_aug = zeros(size(C_aug, 1), size(B_aug, 2));

sys_aug = ss(A_aug, B_aug, C_aug, D_aug);

zhat_ic = zeros(5, 1);
z_aug_ic = [z_ic; zhat_ic];

%% Simulation Setup
t = 0:60:86400;  % 1-minute steps for better numerical stability

% Build external inputs
d_vec = arrayfun(d_t, t);
xa_vec = arrayfun(x_a_dev_t, t);

u_ext = [repmat(d_vec', 1, 5), xa_vec'];

fprintf('Running simulation...\n');
[y_out, t_out, z_aug_traj] = lsim(sys_aug, u_ext, t, z_aug_ic);

%% Extract results
n = size(A, 1);
z_traj = z_aug_traj(:, 1:n);
z_hat_traj = z_aug_traj(:, n+1:end);

% Compute control inputs
u_ctrl_traj = -(z_hat_traj * K');

% Convert to Celsius for plotting
x_traj = z_traj + 20;
x_hat_traj = z_hat_traj + 20;

%% Display results
fprintf('\n=== Performance Summary ===\n');
fprintf('Control effort range: %.1f to %.1f W\n', min(u_ctrl_traj(:)), max(u_ctrl_traj(:)));

if max(abs(u_ctrl_traj(:))) < 1
    fprintf('WARNING: Control effort is virtually zero - system may be unstable\n');
end

%% Plotting
figure('Position', [100, 100, 1200, 1000]);
t_hours = t_out / 3600;

% Plot 1: True room temperatures
subplot(5, 1, 1);
plot(t_hours, x_traj, 'LineWidth', 1.5);
ylabel('Temperature [°C]');
title('True Room Temperatures');
legend('Room 1', 'Room 2', 'Room 3', 'Room 4', 'Room 5');
grid on;
yline(20, 'r--', 'Setpoint');

% Plot 2: Estimated temperatures
subplot(5, 1, 2);
plot(t_hours, x_hat_traj, 'LineWidth', 1.5);
ylabel('Temperature [°C]');
title('Estimated Room Temperatures');
legend('Room 1', 'Room 2', 'Room 3', 'Room 4', 'Room 5');
grid on;
yline(20, 'r--', 'Setpoint');

% Plot 3: Ambient temperature
subplot(5, 1, 3);
ambient_temp = 20 + xa_vec;
plot(t_hours, ambient_temp, 'k', 'LineWidth', 1.5);
ylabel('Temperature [°C]');
title('Ambient Temperature');
legend('Ambient');
grid on;

% Plot 4: Disturbance
subplot(5, 1, 4);
plot(t_hours, d_vec, 'r', 'LineWidth', 1.5);
ylabel('Power [W]');
title('Disturbance Input');
legend('Disturbance');
grid on;

% Plot 5: Control inputs
subplot(5, 1, 5);
plot(t_hours, u_ctrl_traj, 'LineWidth', 1.5);
ylabel('Power [W]');
xlabel('Time [hours]');
title('Control Inputs (Heating)');
legend('Room 1', 'Room 2', 'Room 4');
grid on;

% Check stability
closed_loop_poles = eig(A - B_u*K);
fprintf('Max real part of closed-loop poles: %.6f\n', max(real(closed_loop_poles)));