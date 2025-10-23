nRooms = 5;
nControls = 3;
r = 10*4;
q = 10*2;
x_bar = 20 * ones(nRooms, 1);
T = 86400;  % Fixed: 86400, not 86450
t = linspace(0, T, 1000);  % More points for smoother plots
c = 1/1.08*10^6;

A = [-30 10 0 0 10;
    10 -40 10 0 10;
    0 10 -40 10 10;
    0 0 10 -30 10;
    10 10 10 10 -50];

% FIXED B matrix
B = [1 0 0;
     0 1 0;
     0 0 0;
     0 0 1;
     0 0 0]*c;

C = [0 1 0 0 0;
    0 0 0 1 0];

D = zeros(2, 3);

Q = q * eye(nRooms);
R = r * eye(nControls);
b_ambient = [10; 10; 10; 10; 10];
D_disturbance = [1 0 0 0 0; 
    0 1 0 0 0; 
    0 0 1 0 0; 
    0 0 0 1 0; 
    0 0 0 0 1]*c;

lqrsys = ss(A, B_u, C, D);

[K, P, E] = lqr(lqrsys, Q, R);

x_a = @(t) 20 + 5*sin(2*pi*t/86400);
%% Observer Design (following assignment observability check)
if rank(obsv(A, C)) == size(A,1)
    fprintf('System is observable.\n');
else
    error('System not observable - check C matrix');
end
desired_obs_poles = [-10, -12, -15, -20, -25];
L = place(A', C', desired_obs_poles)';
d_func = @(t) max(0, 500 * sin(2 * pi*t/86400) - 300);
ssign = @(x) (x>0) - (x<0);

% Initial conditions
x0 = zeros(nRooms, 1);
for i = 1:nRooms
    x0(i) = 15 - i * ssign((-1)^i);
end

z0 = x0 - x_bar;
Acl = (A - B*K);

% FIXED: Use z0 for simulation
[t_sim, z_sim] = ode45(@(t,z) buildingDynamics(t, z, A, B, K, x_bar, b_ambient, D_disturbance, x_a, d_func), ...
                       t, z0);

% FIXED: Convert back to actual temperatures
x_sim = z_sim' + x_bar;
x_sim = x_sim';

% FIXED: Control inputs using z_sim
u_sim = zeros(length(t_sim), nControls);

% Then use u_sim_kW in your control plots instead of u_sim
for i = 1:length(t_sim)
    z = z_sim(i,:)';
    u_sim(i,:) = (-K * z)';
end
u_sim_kW = u_sim * 1000;  % Convert Watts to kiloWatts

% ... rest of your plotting code remains the same
legend('Room 1', 'Room 2', 'Room 3', 'Room 4', 'Room 5');
grid on;

% Convert time to hours for plotting
t_hours = t_sim / 3600;

% Create subplots
figure('Position', [100, 100, 1200, 800]);

% Plot 1: Disturbance and ambient temperature
subplot(3,2,1);
d_values = arrayfun(d_func, t_sim);
x_a_values = arrayfun(x_a, t_sim);
plot(t_hours, d_values, 'r-', 'LineWidth', 2); hold on;
plot(t_hours, x_a_values, 'b-', 'LineWidth', 2);
ylabel('Power (W) / Temp (°C)');
title('Disturbance and Ambient Temperature');
legend('d(t) [W]', 'x_a(t) [°C]', 'Location', 'best');
grid on;

% Plot 2: Room temperatures
subplot(3,2,2);
plot(t_hours, x_sim, 'LineWidth', 1.5);
hold on;
yline(20, '--k', 'Setpoint 20°C', 'LineWidth', 1);
ylabel('Temperature (°C)');
title('Room Temperatures');
legend('Room 1', 'Room 2', 'Room 3', 'Room 4', 'Room 5', 'Location', 'best');
grid on;

% Plot 3: Control inputs (in kW)
subplot(3,2,3);
plot(t_hours, u_sim_kW, 'LineWidth', 1.5);
ylabel('Control Input (kW)');
title('Heating Power');
legend('u_1', 'u_2', 'u_3', 'Location', 'best');
grid on;

% Plot 4: Temperature errors
subplot(3,2,4);
errors = x_sim - x_bar';
plot(t_hours, errors, 'LineWidth', 1.5);
ylabel('Temperature Error (°C)');
title('Tracking Errors');
legend('Room 1', 'Room 2', 'Room 3', 'Room 4', 'Room 5', 'Location', 'best');
grid on;

% Plot 5: Control effort over time (in kW)
subplot(3,2,5);
total_power = sum(u_sim_kW, 2);
plot(t_hours, total_power, 'k-', 'LineWidth', 2);
ylabel('Total Power (kW)');
xlabel('Time (hours)');
title('Total Control Effort');
grid on;

% Plot 6: RMS error
subplot(3,2,6);
rms_error = sqrt(mean(errors.^2, 2));
plot(t_hours, rms_error, 'm-', 'LineWidth', 2);
ylabel('RMS Error (°C)');
xlabel('Time (hours)');
title('Overall Tracking Performance');
grid on;

% Debug: Check what's happening
fprintf('Final z values: %.2f, %.2f, %.2f, %.2f, %.2f\n', z_sim(end,:));
fprintf('Final x values: %.2f, %.2f, %.2f, %.2f, %.2f°C\n', x_sim(end,:));
fprintf('Final control inputs: %.2f, %.2f, %.2f W\n', u_sim(end,:));
fprintf('B matrix scaling: %.2e\n', B(1,1));

fprintf('LQR gain matrix K:\n');
disp(K);
fprintf('Max control gain: %.2e\n', max(abs(K(:))));