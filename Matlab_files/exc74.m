%% TEK4090 - Exercise 7.4 - Attitude EKF with varying noise parameters
clear; close all; clc;
rng(42);
q_true = [0; 0; 0; 1]; % Initial attitude quaternion
beta_true = [0.1; 0.1; 0.1] * pi/180/3600; % 0.1 deg/hr -> rad/s
omega_true = [0; 0; pi/2700]; % Rotation rate (90 min orbit)

sigma_star = 0.1 *pi/180; % STAR CAMERA NOISE
sigma_u = sqrt(10) * 1e-10; % GYRO BIAS DRIFT
sigma_v = sqrt(10) * 1e-7;

r1 = [1; 0; 0];
r2 = [0; 1; 0];
r3 = [0; 0; 1];
r = [r1; r2; r3];
% Simulation parameters
dt = 1; % 1 Hz sampling
T = 540; % Reduced for faster testing (9 minutes)
N = T/dt;

% PROBLEM: HOW TO FIGURE OUT WHERE MY SPACESHIP IS POINTING, WHEN MY
% SENSORS ARE NOISY?

function S = skew_symmetric(v)
    S = [0, -v(3), v(2);
         v(3), 0, -v(1);
         -v(2), v(1), 0];
end

function q_out = quat_multiply(q1, q2)
    v1 = q1(1:3); s1 = q1(4);
    v2 = q2(1:3); s2 = q2(4);
    
    q_out = [s1*v2 + s2*v1 + cross(v1, v2);
             s1*s2 - dot(v1, v2)];
    q_out = q_out / norm(q_out);
end

function R = quat_to_dcm(q)
    q0 = q(4); q1 = q(1); q2 = q(2); q3 = q(3);
    R = [1-2*(q2^2+q3^2), 2*(q1*q2-q0*q3), 2*(q1*q3+q0*q2);
         2*(q1*q2+q0*q3), 1-2*(q1^2+q3^2), 2*(q2*q3-q0*q1);
         2*(q1*q3-q0*q2), 2*(q2*q3+q0*q1), 1-2*(q1^2+q2^2)];
end

function q_conj = quat_conjugate(q)
    q_conj = [-q(1:3); q(4)];
end

function q_new = propagate_quartenion(q, omega, dt)
    % Implement discrete version of dq/dt = 0.5 * omega * q
    % Exact solution for constant ω over Δt
    omega_norm = norm(omega);
    
    if omega_norm > 1e-12
        % Axis-angle representation
        axis = omega / omega_norm;
        angle = omega_norm * dt;
        
        % Quaternion for the rotation over Δt
        dq = [axis * sin(angle/2); cos(angle/2)];
        
        % Multiply: q_new = q ⊗ dq
        q_new = quaternion_multiply(q, dq);
    else
        % No rotation
        q_new = q;
    end
    
    % Always normalize
    q_new = q_new / norm(q_new);
end




function q_out = exact_attitude(q0, omega, t)
    omega_norm = norm(omega);
    if omega_norm > 1e-12
        axis = omega / omega_norm;
        angle = omega_norm * t;
        q_out = [axis * sin(angle/2); cos(angle/2)];
        q_out = quat_multiply(q0, q_out);
    else
        q_out = q0;
    end
    q_out = q_out / norm(q_out);
end


% True system: ω_true = constant
% Measurements: ω_meas = ω_true + β + η_v
function gyro_meas = generate_gyro_measurements(omega_true, beta_true, sigma_v, N)
    gyro_meas = zeros(3, N);
    for k = 1:N
        eta_v = sigma_v * randn(3,1);  % Measurement noise
        gyro_meas(:,k) = omega_true + beta_true(:,k) + eta_v;
    end
end

function star_meas = generate_star_measurements(q_true, sigma_star, N)
    % Reference vectors in inertial frame
    r1 = [1;0;0]; r2 = [0;1;0]; r3 = [0;0;1];
    
    star_meas = zeros(9, N);
    for k = 1:N
        R = quat_to_dcm(q_true(:,k));  % True rotation matrix
        noise = sigma_star * randn(9,1);
        
        % Measured body vectors
        b1 = R * r1 + noise(1:3);
        b2 = R * r2 + noise(4:6); 
        b3 = R * r3 + noise(7:9);
        
        star_meas(:,k) = [b1/norm(b1); b2/norm(b2); b3/norm(b3)];
    end
end

function [q_est, beta_est, P] = attitude_ekf(q0, beta0, P0, gyro_meas, star_meas, dt, sigma_v, sigma_u, sigma_star)
    % Initialize
    N = size(gyro_meas, 2);
    q_est = zeros(4, N);
    beta_est = zeros(3, N);
    P = zeros(6, 6, N);
    
    q_est(:,1) = q0;
    beta_est(:,1) = beta0;
    P(:,:,1) = P0;
    
    % Reference vectors for star measurements
    r1 = [1;0;0]; r2 = [0;1;0]; r3 = [0;0;1];
    r_vec = [r1; r2; r3];
    
    for k = 1:N-1
        %% Propagation Step
        % Estimated angular velocity
        omega_est = gyro_meas(:,k) - beta_est(:,k);
        
        % Quaternion propagation (discrete)
        omega_norm = norm(omega_est);
        if omega_norm > 1e-12
            psi = sin(0.5*omega_norm*dt) * omega_est / omega_norm;
            cos_term = cos(0.5*omega_norm*dt);            
            Omega = [cos_term*eye(3) - skew_symmetric(psi),  psi;
                     -psi',                              cos_term];
            q_est(:,k+1) = Omega * q_est(:,k);
        else
            q_est(:,k+1) = q_est(:,k);
        end
        q_est(:,k+1) = q_est(:,k+1) / norm(q_est(:,k+1));
        
        % Bias propagation (constant)
        beta_est(:,k+1) = beta_est(:,k);
        
        % Covariance propagation
        F = [-skew_symmetric(omega_est), -eye(3);
              zeros(3,3),               zeros(3,3)];
        
        % Discrete state transition matrix (simplified - you can implement Eq. 7.45)
        Phi = eye(6) + F * dt;
        
        G = [-eye(3), zeros(3,3);
              zeros(3,3), eye(3,3)];
        
        Qc = [sigma_v^2 * eye(3), zeros(3,3);
              zeros(3,3), sigma_u^2 * eye(3)];
        Qd = G * Qc * G' * dt;
        
        P_pred = Phi * P(:,:,k) * Phi' + Qd;
        
        %% Update Step (when star measurements available)
        if ~isempty(star_meas)
            % Current propagated estimates
            q_pred = q_est(:,k+1);
            R_pred = quat_to_dcm(q_pred);
            
            % Build sensitivity matrix H
            H = zeros(9, 6);
            H(1:3, 1:3) = skew_symmetric(R_pred * r1);
            H(4:6, 1:3) = skew_symmetric(R_pred * r2);
            H(7:9, 1:3) = skew_symmetric(R_pred * r3);
            % Last 3 columns remain zero (no direct bias measurement)
            
            % Measurement covariance
            R_meas = sigma_star^2 * eye(9);
            
            % Kalman gain
            K = P_pred * H' * inv(H * P_pred * H' + R_meas);
            
            % Predicted measurements
            b1_pred = R_pred * r1;
            b2_pred = R_pred * r2;
            b3_pred = R_pred * r3;
            y_pred = [b1_pred; b2_pred; b3_pred];
            
            % Measurement residual
            residual = star_meas(:,k+1) - y_pred;
            
            % Error state update
            delta_x = K * residual;
            delta_alpha = delta_x(1:3);
            delta_beta = delta_x(4:6);
            
            % Update covariance
            P_upd = (eye(6) - K * H) * P_pred;
            P(:,:,k+1) = P_upd;
            
            % Update quaternion (multiplicative)
            delta_q = [0.5 * delta_alpha; 1];
            delta_q = delta_q / norm(delta_q);
            q_est(:,k+1) = quat_multiply(delta_q, q_pred);
            q_est(:,k+1) = q_est(:,k+1) / norm(q_est(:,k+1));
            
            % Update bias
            beta_est(:,k+1) = beta_est(:,k+1) + delta_beta;
        else
            P(:,:,k+1) = P_pred;
        end
    end
end

%% Full simulation with different noise scenarios
scenarios = {
    struct('sigma_u', sqrt(10)*1e-10, 'sigma_v', sqrt(10)*1e-7, 'name', 'Baseline'),
    struct('sigma_u', sqrt(10)*1e-10, 'sigma_v', sqrt(10)*1e-5, 'name', 'High_sigma_v'), 
    struct('sigma_u', sqrt(10)*1e-8,  'sigma_v', sqrt(10)*1e-7, 'name', 'High_sigma_u')
};

% Simulation parameters
dt = 1; % 1 Hz sampling  
T = 5400; % 90 minutes in seconds (full orbit)
N = T/dt;

% True angular velocity (90 min orbit = 0.00116 rad/s)
omega_true = [0; 0; 2*pi/(90*60)]; % Rotation about z-axis

% Extended simulation
results = cell(1,3);

for scenario_idx = 1:3
    scenario = scenarios{scenario_idx};
    fprintf('\n=== Running scenario: %s ===\n', scenario.name);
    
    % Generate true trajectory
    q_true_history = zeros(4, N);
    beta_true_history = zeros(3, N);
    
    q_true_history(:,1) = [0; 0; 0; 1];
    beta_true_history(:,1) = [0.1; 0.1; 0.1] * pi/180/3600;
    
    % Constant bias for true system
    for k = 2:N
        q_true_history(:,k) = exact_attitude(q_true_history(:,k-1), omega_true, dt);
        beta_true_history(:,k) = beta_true_history(:,1); % Constant true bias
    end
    
    % Generate measurements
    gyro_meas = zeros(3, N);
    for k = 1:N
        eta_v = scenario.sigma_v * randn(3,1);
        gyro_meas(:,k) = omega_true + beta_true_history(:,k) + eta_v;
    end
    
    star_meas = zeros(9, N);
    for k = 1:N
        R_true = quat_to_dcm(q_true_history(:,k));
        noise = sigma_star * randn(9,1);
        
        b1 = R_true * r1 + noise(1:3);
        b2 = R_true * r2 + noise(4:6);
        b3 = R_true * r3 + noise(7:9);
        
        star_meas(:,k) = [b1/norm(b1); b2/norm(b2); b3/norm(b3)];
    end
    
    % Run EKF
    q0_est = q_true_history(:,1); % Start with perfect knowledge
    beta0_est = [0; 0; 0]; % Zero initial bias estimate
    
    P0a = 3.0462e-6;
    P0b = 9.4018e-13;
    P0 = diag([P0a, P0a, P0a, P0b, P0b, P0b]);
    
    [q_est, beta_est, P] = attitude_ekf(q0_est, beta0_est, P0, gyro_meas, star_meas, dt, ...
                                       scenario.sigma_v, scenario.sigma_u, sigma_star);
    
    % Compute errors
    attitude_errors = zeros(3, N);
    bias_errors = zeros(3, N);
    
    for k = 1:N
        % Attitude error
        q_err = quat_multiply(quat_conjugate(q_est(:,k)), q_true_history(:,k));
        angle_err = 2 * acos(abs(q_err(4))) * 180/pi; % degrees
        attitude_errors(:,k) = [angle_err; 0; 0]; % Store scalar error
        
        % Bias error  
        bias_errors(:,k) = (beta_est(:,k) - beta_true_history(:,k)) * (180/pi * 3600); % deg/hr
    end
    
    % Store results
    results{scenario_idx} = struct();
    results{scenario_idx}.name = scenario.name;
    results{scenario_idx}.q_est = q_est;
    results{scenario_idx}.beta_est = beta_est;
    results{scenario_idx}.attitude_errors = attitude_errors;
    results{scenario_idx}.bias_errors = bias_errors;
    results{scenario_idx}.P = P;
    
    fprintf('Final attitude error: %.4f deg\n', attitude_errors(1,end));
    fprintf('Final bias errors: [%.3f, %.3f, %.3f] deg/hr\n', bias_errors(:,end));
end

%% Plot results like Figure 7.2 and 7.3
time_minutes = (0:N-1)*dt/60;

figure('Position', [100, 100, 1200, 800]);

% Attitude errors
subplot(2,2,1);
hold on;
colors = ['b', 'r', 'g'];
for i = 1:3
    plot(time_minutes, results{i}.attitude_errors(1,:)*60, colors(i), 'LineWidth', 1.5);
end
ylabel('Attitude Error (arcmin)');
xlabel('Time (min)');
title('Attitude Estimation Errors');
legend({scenarios{1}.name, scenarios{2}.name, scenarios{3}.name}, 'Location', 'best');
grid on;

% Bias estimates
subplot(2,2,2);
hold on;
for i = 1:3
    bias_deg_hr = results{i}.beta_est * (180/pi * 3600);
    plot(time_minutes, bias_deg_hr(1,:), colors(i), 'LineWidth', 1.5);
end
plot([0, time_minutes(end)], [0.1, 0.1], 'k--', 'LineWidth', 2, 'DisplayName', 'True Bias');
ylabel('X-Bias Estimate (deg/hr)');
xlabel('Time (min)');
title('Gyro Bias Estimation');
legend('Location', 'best');
grid on;

% 3σ bounds for baseline scenario
subplot(2,2,3);
attitude_std = sqrt(squeeze(P(1,1,:))) * 180/pi; % Convert to degrees
three_sigma = 3 * attitude_std * 60; % Convert to arcmin

plot(time_minutes, results{1}.attitude_errors(1,:)*60, 'b', 'LineWidth', 1.5);
hold on;
plot(time_minutes, three_sigma, 'r--', 'LineWidth', 1.5);
plot(time_minutes, -three_sigma, 'r--', 'LineWidth', 1.5);
ylabel('Attitude Error (arcmin)');
xlabel('Time (min)');
title('Attitude Errors with 3σ Bounds (Baseline)');
legend('Attitude Error', '3σ Bound', 'Location', 'best');
grid on;

% All bias components for baseline
subplot(2,2,4);
bias_deg_hr = results{1}.beta_est * (180/pi * 3600);
plot(time_minutes, bias_deg_hr(1,:), 'b', 'LineWidth', 1.5);
hold on;
plot(time_minutes, bias_deg_hr(2,:), 'r', 'LineWidth', 1.5);
plot(time_minutes, bias_deg_hr(3,:), 'g', 'LineWidth', 1.5);
plot([0, time_minutes(end)], [0.1, 0.1], 'k--', 'LineWidth', 2);
ylabel('Bias Estimates (deg/hr)');
xlabel('Time (min)');
title('All Bias Components (Baseline)');
legend('X-bias', 'Y-bias', 'Z-bias', 'True', 'Location', 'best');
grid on;

fprintf('\n=== Summary ===\n');
for i = 1:3
    final_att_err = results{i}.attitude_errors(1,end);
    final_bias_err = mean(abs(results{i}.bias_errors(:,end)));
    fprintf('%s: Final attitude error = %.3f arcmin, Final bias error = %.3f deg/hr\n', ...
            scenarios{i}.name, final_att_err*60, final_bias_err);
end