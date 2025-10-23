%% TEK4090 - Assignment 2, Exercise 4: Building Temperature Control
% Implementation following the approach from the assignment solution

clear all; close all; clc;

%% Define model parameters (following assignment notation)
c = 1/(1.08*10^6);  % [K/J] - Inverse thermal capacitance

% Setpoint temperature 
x_bar = [0; 0; 0; 0; 0] + 20 + 273;  % [K] - Setpoint in Kelvin

% Initial conditions (converted to Kelvin as in assignment)
x_ic = [16; 13; 18; 11; 20] + 273;  % [K]

% z = x - x_bar (deviation coordinates as in assignment)
z_ic = x_ic - x_bar;  % [K]

function d_t = dtr(t)
    d = 500*sin(2*pi*t/86400) - 300;
    d_t = max(0, d);
end
function xa_t = amb(t)
    xa_t = 5*sin(2*pi*t/86400) - 300;
end

%% System matrices (following assignment structure)
A = [-30,  10,   0,   0,  10;
      10, -40,  10,   0,  10;
       0,  10, -40,  10,  10;
       0,   0,  10, -30,  10;
      10,  10,  10,  10, -50];

B_u = [c, 0, 0;
       0, c, 0; 
       0, 0, 0;
       0, 0, c;
       0, 0, 0];

B_d = eye(5) * c;
B_a = [10; 10; 10; 10; 10] * c;  % [1/s]
C = [1, 0, 0, 0, 0;
     0, 0, 1, 0, 0];
D = zeros(2, 3);

%% LQR Controller Design (using assignment approach)
Q = eye(5) * 1e20;    % State penalty
R = eye(3) * 1e-10;   % Control penalty
lqr_sys = ss(A, B_u, C, D);
[K, S, P] = lqr(lqr_sys, Q, R);

%% Observer Design (following assignment observability check)
if rank(obsv(A, C)) == size(A,1)
    fprintf('System is observable.\n');
else
    error('System not observable - check C matrix');
end
desired_obs_poles = [-10, -12, -15, -20, -25];
L = place(A', C', desired_obs_poles)';

%% Augmented System for Simulation (following assignment structure)
A_aug = [A,          -B_u*K;
         L*C,   A - B_u*K - L*C];

B_ext_top = [B_d, B_a];
B_ext_bottom = [B_d, B_a];
B_aug = [B_ext_top; B_ext_bottom];

C_aug = [C, zeros(size(C))];
D_aug = zeros(size(C_aug, 1), size(B_aug, 2));

sys_aug = ss(A_aug, B_aug, C_aug, D_aug);

zhat_ic = zeros(size(z_ic));  % Initial estimate
z_aug_ic = [z_ic; zhat_ic];

%% Simulation Setup (following assignment time handling)
t = 0:1:86400;  % 10-second steps for efficiency

d_vec = arrayfun(@dtr, t(:));
xa_vec = arrayfun(@amb, t(:));

u_ext = [repmat(d_vec, 1, 5), xa_vec];
[y_out, t_out, z_aug_traj] = lsim(sys_aug, u_ext, t, z_aug_ic);

n = size(A, 1);
z_traj = z_aug_traj(:, 1:n);           % True states
z_hat_traj = z_aug_traj(:, n+1:end);   % Estimated states

y_traj = y_out;

% Compute control inputs u(t) = -K * z_hat(t)
u_ctrl_traj = -(z_hat_traj * K');

% Convert back to actual temperatures in Celsius
%x_traj = z_traj + (x_bar - 273);      % Convert to °C
%x_hat_traj = z_hat_traj + (x_bar - 273); % Convert to °C



%% Comprehensive Plotting (improved version of assignment plots)
figure;

% Plot 1: True room temperatures (°C)
subplot(5, 1, 1);
plot(t_out, z_traj + 20, 'LineWidth', 1.5);
ylabel('Temp [°C]');
title('True Room Temperatures');
legend('Room 1', 'Room 2', 'Room 3', 'Room 4', 'Room 5', 'Location', 'best');
grid on;
yline(20, 'r--', 'Setpoint', 'LineWidth', 1);

% Plot 2: Estimated temperatures (°C)
subplot(5, 1, 2);
plot(t_out, z_hat_traj + 20, 'LineWidth', 1.5);
ylabel('Temp [°C]');
title('Estimated Room Temperatures');
legend('Room 1', 'Room 2', 'Room 3', 'Room 4', 'Room 5', 'Location', 'best');
grid on;
yline(20, 'r--', 'Setpoint', 'LineWidth', 1);

% Plot 3: Ambient temperature (°C)
subplot(5, 1, 3);
plot(t_out, xa_vec + 20, 'k', 'LineWidth', 1.5); ylabel('Temp [°C]'); title('Ambient Temperature');
legend('xa'); grid on;

% Plot 4: Disturbance profile
subplot(5, 1, 4);
plot(t_out, repmat(d_vec, 1, 5), 'r', 'LineWidth', 1.5); ylabel('Power [W]'); title('Disturbance Input'); 
legend('d1', 'd2', 'd3', 'd4', 'd5'); grid on;

% Plot 5: Control inputs
subplot(5, 1, 5);
plot(t_out, u_ctrl_traj, 'LineWidth', 1.5); ylabel('Power [W]'); xlabel('Time [hours]');
title('Control Inputs (Heating)'); legend('u1', 'u2', 'u3'); grid on;


% Check if LQR produced stable closed-loop poles
closed_loop_poles = eig(A - B_u*K);
fprintf('Max real part of closed-loop poles: %.3f\n', max(real(closed_loop_poles)));
if any(real(closed_loop_poles) >= 0)
    error('LQR produced unstable closed-loop system!');
end