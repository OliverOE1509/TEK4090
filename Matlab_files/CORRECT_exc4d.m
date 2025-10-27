nRooms = 5;
nControls = 3;
nMeasurements = 2;

q = 1e4;  
r = 1e2;  
x_bar = 20 * ones(nRooms, 1);
sigma_meas = 0.1;
T = 86400;
t = linspace(0, T, 1000);
c = 1/(1.08e6);

A = [-30 10 0 0 10;
     10 -40 10 0 10;
     0 10 -40 10 10;
     0 0 10 -30 10;
     10 10 10 10 -50] ;

B = [1 0 0;
     0 1 0;
     0 0 0;
     0 0 1;
     0 0 0] * c;

C = [0 1 0 0 0;   
     0 0 0 1 0];   

Q = q * eye(nRooms);
R = r * eye(nControls);
[K, P, E] = lqr(A, B, Q, R);

if rank(obsv(A, C)) == nRooms
    fprintf('System is observable.\n');
else
    fprintf('System is NOT observable\n');
end

obs_poles_multiplier = 2;  
desired_obs_poles = obs_poles_multiplier * max(real(E)) * linspace(1, 1.5, nRooms);
L = place(A', C', desired_obs_poles)';

B_ambient = [10; 10; 10; 10; 10] * c;
D_disturbance = eye(5) * c;
x_a = @(t) 20 + 5*sin(2*pi*t/86400);
d_func = @(t) max(0, 500 * sin(2 * pi*t/86400) - 300);

ssign = @(x) (x>0) - (x<0);
x0_true = zeros(nRooms, 1);
for i = 1:nRooms
    x0_true(i) = 15 - i * ssign((-1)^i);
end
z0_true = x0_true - x_bar;

x0_est = 18 * ones(nRooms, 1);  
z0_est = x0_est - x_bar;
disp(z0_est)

z0_aug = [z0_true; z0_est];

function dz_aug = augmentedDynamics(t, z_aug, A, B, K, L, C, B_ambient, D_disturbance, x_a, d_func, sigma_meas)
    z_true = z_aug(1:5);
    z_est = z_aug(6:10);
    
    y_meas = C * z_true + sigma_meas * randn(2,1);
    
    dz_true = A * z_true + B * (-K * z_est) + ...
              B_ambient * (x_a(t) - 20) + ...
              D_disturbance * d_func(t) * ones(5,1);
    
    dz_est = A * z_est + B * (-K * z_est) + L * (y_meas - C * z_est);
    
    dz_aug = [dz_true; dz_est];
end

[t_sim, z_aug] = ode45(@(t,z) augmentedDynamics(t, z, A, B, K, L, C, B_ambient, D_disturbance, x_a, d_func, sigma_meas), ...
                       t, z0_aug);
z_true = z_aug(:, 1:5)';
z_est = z_aug(:, 6:10)';

x_true = z_true + x_bar;
x_est = z_est + x_bar;

u_sim = zeros(length(t_sim), nControls);
for i = 1:length(t_sim)
    u_sim(i,:) = (-K * z_est(:,i))';
end


t_hours = t_sim / 3600;

figure('Position', [100, 100, 1400, 1000]);

subplot(3, 2, 1);
d_values = arrayfun(d_func, t_sim);
x_a_values = arrayfun(x_a, t_sim);
plot(t_hours, d_values, 'r-', 'LineWidth', 2); hold on;
plot(t_hours, x_a_values, 'b-', 'LineWidth', 2);
ylabel('Power (W) / Temp (°C)');
title('External Influences');
legend('Disturbance d(t) [W]', 'Ambient Temp x_a(t) [°C]', 'Location', 'best');
grid on;

subplot(3, 2, 2);
plot(t_hours, x_true', 'LineWidth', 2); hold on;
plot(t_hours, x_est', '--', 'LineWidth', 1);
yline(20, '--k', 'Setpoint 20°C', 'LineWidth', 1);
ylabel('Temperature (°C)');
title('True (solid) vs Estimated (dashed) Temperatures');
legend('Room 1', 'Room 2', 'Room 3', 'Room 4', 'Room 5', 'Location', 'best');
grid on;

subplot(3, 2, 3);
plot(t_hours, u_sim, 'LineWidth', 2);
ylabel('Heating Power (W)');
title('Control Inputs');
legend('Room 1 Heater', 'Room 2 Heater', 'Room 4 Heater', 'Location', 'best');
grid on;

subplot(3, 2, 4);
est_errors = x_true' - x_est';
plot(t_hours, est_errors, 'LineWidth', 1.5);
ylabel('Estimation Error (°C)');
title('Observer Performance: True - Estimated');
legend('Room 1', 'Room 2', 'Room 3', 'Room 4', 'Room 5', 'Location', 'best');
grid on;

subplot(3, 2, 5);
total_power = sum(u_sim, 2);
plot(t_hours, total_power, 'k-', 'LineWidth', 2);
ylabel('Total Power (W)');
xlabel('Time (hours)');
title('Total Control Effort');
grid on;

subplot(3, 2, 6);
tracking_errors = x_true' - x_bar';
plot(t_hours, tracking_errors, 'LineWidth', 1.5);
ylabel('Tracking Error (°C)');
xlabel('Time (hours)');
title('True System Tracking Performance');
legend('Room 1', 'Room 2', 'Room 3', 'Room 4', 'Room 5', 'Location', 'best');
grid on;
