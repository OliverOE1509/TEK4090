clear all; close all; clc;

% System parameters - FIXED SCALING
C = 1.08e6;
A = [-30 10 0 0 10;
     10 -40 10 0 10;
     0 10 -40 10 10;
     0 0 10 -30 10;
     10 10 10 10 -50] / C;  % SCALE THE A MATRIX by thermal capacitance

B_u = [1/C 0 0;
       0 1/C 0;
       0 0 0;
       0 0 1/C;
       0 0 0];

B_a = [10;10;10;10;10];
B_d = eye(5,1) / C;

% LQR weights - ADJUSTED FOR BETTER PERFORMANCE
Q = eye(5) * 1e6;  % Increased to prioritize temperature tracking
R = eye(3) * 1e-3; % Decreased to allow more control action
B = [B_u, B_d, B_a];

C_output = [1 0 0 0 0;
    0 1 0 0 0;];

D = zeros(2,3);

LQR_sys = ss(A, B_u, C_output, D);
% Compute LQR gain
[K, S, E] = lqr(LQR_sys, Q, R);

%% Simulation with disturbances and ambient temperature
% Initial conditions
x0 = [18; 17; 16; 19; 15];  % Different initial temperatures
x_bar = 20;                  % Setpoint temperature
z0 = x0 - x_bar;             % Deviation from setpoint

% Time vector (24 hours in seconds, 1-minute steps)
t = 0:60:86400;  

% Disturbance: internal heat gains (people, equipment)
d_vec = max(0, 300 * sin(2 * pi * t / 86400) - 150);  

% Ambient temperature: daily variation
x_a_full = 20 + 8 * sin(2 * pi * (t - 6*3600) / 86400);  % Peak at 2 PM

% Preallocate
z_traj = zeros(length(t), 5);
u_traj = zeros(length(t), 3);
x_traj = zeros(length(t), 5);

z_current = z0;

% Simulation loop
for k = 1:length(t)
    z_traj(k,:) = z_current';
    x_traj(k,:) = z_current' + x_bar;
    
    % Control law: u = -K * z
    u_current = -K * z_current;
    u_traj(k,:) = u_current';
    
    % System dynamics
    if k < length(t)
        disturbance_effect = B_d * d_vec(k);
        ambient_effect = B_a * (x_a_full(k) - 20) / C;
        
        % Euler integration (60s step)
        dz = A * z_current + B_u * u_current + disturbance_effect + ambient_effect;
        z_current = z_current + dz * 60;
    end
end

%% Plot results
figure('Position', [100, 100, 1200, 800]);

% Plot 1: Room temperatures
subplot(3,1,1);
plot(t/3600, x_traj, 'LineWidth', 2);
ylabel('Temperature [°C]');
title('Room Temperatures with LQR Control');
legend('Room 1', 'Room 2', 'Room 3', 'Room 4', 'Room 5', 'Location', 'best');
grid on;
yline(20, 'r--', 'Setpoint 20°C', 'LineWidth', 1.5);
xlim([0 24]);
ylim([18 22]);

% Plot 2: Control inputs
subplot(3,1,2);
plot(t/3600, u_traj/1000, 'LineWidth', 2);
ylabel('Heating Power [kW]');
title('Control Inputs');
legend('Room 1', 'Room 2', 'Room 4', 'Location', 'best');
grid on;
xlim([0 24]);

% Plot 3: Disturbances
subplot(3,1,3);
yyaxis left
plot(t/3600, d_vec, 'LineWidth', 2);
ylabel('Internal Heat [W]');
yyaxis right
plot(t/3600, x_a_full, 'LineWidth', 2);
ylabel('Ambient Temp [°C]');
xlabel('Time [hours]');
title('Disturbances');
legend('Internal Heat Gain', 'Ambient Temperature', 'Location', 'best');
grid on;
xlim([0 24]);

%% Display results
fprintf('\n=== LQR CONTROL SYSTEM PERFORMANCE ===\n');
fprintf('Initial temperatures:  [%.1f, %.1f, %.1f, %.1f, %.1f]°C\n', x0);
fprintf('Final temperatures:    [%.2f, %.2f, %.2f, %.2f, %.2f]°C\n', x_traj(end,:));
fprintf('Temperature RMS error: [%.3f, %.3f, %.3f, %.3f, %.3f]°C\n', ...
    rms(x_traj - 20));
fprintf('Maximum control power: [%.1f, %.1f, %.1f] W\n', max(abs(u_traj)));
fprintf('Average control power: [%.1f, %.1f, %.1f] W\n', mean(abs(u_traj)));

% Check controlability and observability
fprintf('\n=== SYSTEM PROPERTIES ===\n');
Co = ctrb(A, B_u);
Ob = obsv(A, C_output);
fprintf('Controlability rank: %d (full rank = %d)\n', rank(Co), size(A,1));
fprintf('Observability rank: %d (full rank = %d)\n', rank(Ob), size(A,1));

% LQR properties
fprintf('LQR poles real parts: ');
fprintf('%.3f ', real(E));
fprintf('\nAll poles stable? %s\n', mat2str(all(real(E) < 0)));

%% Additional analysis: Steady-state performance
steady_state_error = mean(x_traj(end-100:end,:) - 20);
fprintf('\n=== STEADY-STATE ANALYSIS ===\n');
fprintf('Steady-state errors: [%.4f, %.4f, %.4f, %.4f, %.4f]°C\n', steady_state_error);
fprintf('Max steady-state error: %.4f°C\n', max(abs(steady_state_error)));

% Plot steady-state region
figure('Position', [100, 100, 800, 400]);
plot(t/3600, x_traj - 20, 'LineWidth', 1.5);
ylabel('Temperature Error [°C]');
xlabel('Time [hours]');
title('Temperature Deviation from Setpoint');
legend('Room 1', 'Room 2', 'Room 3', 'Room 4', 'Room 5', 'Location', 'best');
grid on;
xlim([20 24]);  % Last 4 hours
yline(0, 'k--', 'Setpoint', 'LineWidth', 2);
ylim([-0.5 0.5]);