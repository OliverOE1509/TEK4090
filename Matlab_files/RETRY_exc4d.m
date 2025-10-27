nRooms = 5;
nControls = 3;
r = 1e2;
q = 1e6;
x_bar = 20 * ones(nRooms, 1);
T = 86400;
t = linspace(0, T, 1000);
c = 1/(1.08e6);  % FIXED: Proper scaling factor

% FIXED: Correct A matrix with proper scaling
A = [-30 10 0 0 10;
     10 -40 10 0 10;
     0 10 -40 10 10;
     0 0 10 -30 10;
     10 10 10 10 -50] * c;  % SCALE A matrix by c

% FIXED: B matrix with proper scaling
B = [1 0 0;
     0 1 0;
     0 0 0;
     0 0 1;
     0 0 0] * c;

C = [0 1 0 0 0;
     0 0 0 1 0];

D = zeros(2, 3);

Q = q * eye(nRooms);
R = r * eye(nControls);

% FIXED: Ambient coupling matrix (scaled by c)
B_ambient = [10; 10; 10; 10; 10] * c;

% FIXED: Disturbance matrix (scaled by c)
D_disturbance = eye(5) * c;

% FIXED: Use continuous-time LQR for the A, B matrices
[K, P, E] = lqr(A, B, Q, R);

% Ambient temperature function
x_a = @(t) 20 + 5*sin(2*pi*t/86400);

% Observer Design
if rank(obsv(A, C)) == size(A,1)
    fprintf('System is observable.\n');
else
    error('System not observable - check C matrix');
end

desired_obs_poles = [-0.1, -0.12, -0.15, -0.2, -0.25];  % FIXED: Slower poles for stability
L = place(A', C', desired_obs_poles)';

% Disturbance function
d_func = @(t) max(0, 500 * sin(2 * pi*t/86400) - 300);

% Sign function
ssign = @(x) (x>0) - (x<0);

% Initial conditions
x0 = zeros(nRooms, 1);
for i = 1:nRooms
    x0(i) = 15 - i * ssign((-1)^i);
end

z0 = x0 - x_bar;

% FIXED: Dynamics function that properly handles all inputs
buildingDynamics = @(t, z) A*z + B*(-K*z) + B_ambient*(x_a(t) - 20) + D_disturbance*d_func(t)*ones(5,1);

% Simulate the system
[t_sim, z_sim] = ode45(@(t,z) buildingDynamics(t, z), t, z0);

% Convert back to actual temperatures
x_sim = z_sim + x_bar';  % FIXED: Proper dimensions

% Calculate control inputs
u_sim = zeros(length(t_sim), nControls);
for i = 1:length(t_sim)
    z = z_sim(i,:)';
    u_sim(i,:) = (-K * z)';
end

% Convert time to hours for plotting
t_hours = t_sim / 3600;

% Create comprehensive plots
figure('Position', [100, 100, 1400, 900]);

% Plot 1: Disturbance and ambient temperature
subplot(3, 2, 1);
d_values = arrayfun(d_func, t_sim);
x_a_values = arrayfun(x_a, t_sim);
plot(t_hours, d_values, 'r-', 'LineWidth', 2); hold on;
plot(t_hours, x_a_values, 'b-', 'LineWidth', 2);
ylabel('Power (W) / Temp (°C)');
title('External Influences');
legend('Disturbance d(t) [W]', 'Ambient Temp x_a(t) [°C]', 'Location', 'best');
grid on;
xlabel('Time (hours)');

% Plot 2: Room temperatures
subplot(3, 2, 2);
plot(t_hours, x_sim, 'LineWidth', 2);
hold on;
yline(20, '--k', 'Setpoint 20°C', 'LineWidth', 1, 'LabelVerticalAlignment', 'bottom');
ylabel('Temperature (°C)');
title('Room Temperatures');
legend('Room 1', 'Room 2', 'Room 3', 'Room 4', 'Room 5', 'Location', 'best');
grid on;
xlabel('Time (hours)');

% Plot 3: Control inputs (in Watts)
subplot(3, 2, 3);
plot(t_hours, u_sim, 'LineWidth', 2);
ylabel('Heating Power (W)');
title('Control Inputs');
legend('Room 1 Heater', 'Room 2 Heater', 'Room 4 Heater', 'Location', 'best');
grid on;
xlabel('Time (hours)');

% Plot 4: Temperature errors
subplot(3, 2, 4);
errors = x_sim - x_bar';
plot(t_hours, errors, 'LineWidth', 1.5);
ylabel('Temperature Error (°C)');
title('Tracking Errors');
legend('Room 1', 'Room 2', 'Room 3', 'Room 4', 'Room 5', 'Location', 'best');
grid on;
xlabel('Time (hours)');

% Plot 5: Total control effort
subplot(3, 2, 5);
total_power = sum(u_sim, 2);
plot(t_hours, total_power, 'k-', 'LineWidth', 2);
ylabel('Total Power (W)');
xlabel('Time (hours)');
title('Total Control Effort');
grid on;

% Plot 6: Control effort distribution
subplot(3, 2, 6);
avg_power = mean(u_sim, 1);
bar(avg_power, 'FaceAlpha', 0.7);
ylabel('Average Power (W)');
title('Average Control Distribution');
xticklabels({'Room 1', 'Room 2', 'Room 4'});
grid on;

% Display key results
fprintf('\n=== SIMULATION RESULTS ===\n');
fprintf('Final temperatures: %.2f, %.2f, %.2f, %.2f, %.2f°C\n', x_sim(end,:));
fprintf('Final errors: %.2f, %.2f, %.2f, %.2f, %.2f°C\n', errors(end,:));
fprintf('Final control inputs: %.1f, %.1f, %.1f W\n', u_sim(end,:));
fprintf('Average control inputs: %.1f, %.1f, %.1f W\n', avg_power);
fprintf('Maximum total power: %.1f W\n', max(total_power));

% Additional analysis: Check if temperatures are reasonable
fprintf('\n=== PERFORMANCE ANALYSIS ===\n');
final_rms_error = sqrt(mean(errors(end,:).^2));
fprintf('Final RMS error: %.3f°C\n', final_rms_error);

% Check if temperatures stay within reasonable bounds
max_temp = max(x_sim(:));
min_temp = min(x_sim(:));
fprintf('Temperature range: %.1f°C to %.1f°C\n', min_temp, max_temp);

% Analyze control activity
total_energy = trapz(t_sim, total_power);
fprintf('Total energy used: %.0f J\n', total_energy);