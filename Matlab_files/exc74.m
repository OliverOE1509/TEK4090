%% TEK4090 - Exercise 7.4 - Attitude EKF with varying noise parameters
clear; close all; clc;

rng(42); % For reproducible results

%% Parameters from Example 7.1
sigma_star = 0.01; % Star camera noise (rad)
sigma_u = sqrt(10)*1e-10; % Gyro bias drift (rad/s^1.5)  
sigma_v = sqrt(10)*1e-7;  % Gyro measurement noise (rad/s^0.5)

% Initial conditions
q_true = [0; 0; 0; 1]; % Initial attitude quaternion
beta_true = [0.1; 0.1; 0.1] * pi/180/3600; % 0.1 deg/hr -> rad/s
omega_true = [0; 0; pi/2700]; % Rotation rate (90 min orbit)

% Simulation parameters
dt = 1; % 1 Hz sampling
T = 540; % Reduced for faster testing (9 minutes)
N = T/dt;

%% Test different noise scenarios
scenarios = {
    struct('name', 'Low Noise', 'sigma_u', sigma_u, 'sigma_v', sigma_v),
    struct('name', 'High sigma_v', 'sigma_u', sigma_u, 'sigma_v', 1e-5),
    struct('name', 'High sigma_u', 'sigma_u', 1e-8, 'sigma_v', sigma_v)
};

results = cell(1, length(scenarios));

for s = 1:length(scenarios)
    fprintf('Running scenario %d: %s\n', s, scenarios{s}.name);
    results{s} = run_attitude_ekf(scenarios{s}, q_true, beta_true, omega_true, dt, N, sigma_star);
end

%% Plot results
figure('Position', [100, 100, 1200, 800]);

% Attitude error
subplot(3,1,1);
hold on;
colors = ['b', 'r', 'g'];
for s = 1:length(scenarios)
    plot(results{s}.time, results{s}.attitude_error, colors(s), 'LineWidth', 1.5);
end
ylabel('Attitude Error (rad)');
title('Attitude Estimation Error');
legend({scenarios{1}.name, scenarios{2}.name, scenarios{3}.name}, 'Location', 'northeast');
grid on;

% Bias estimation error
subplot(3,1,2);
hold on;
for s = 1:length(scenarios)
    plot(results{s}.time, results{s}.bias_error, colors(s), 'LineWidth', 1.5);
end
ylabel('Bias Error (rad/s)');
title('Gyro Bias Estimation Error');
grid on;

% 3-sigma bounds coverage
subplot(3,1,3);
hold on;
for s = 1:length(scenarios)
    plot(results{s}.time, results{s}.coverage, colors(s), 'LineWidth', 1.5);
end
ylabel('Error within 3σ (%)');
xlabel('Time (s)');
title('Estimation Error within 3σ Bounds');
grid on;

%% Performance summary
fprintf('\n=== PERFORMANCE SUMMARY ===\n');
fprintf('Scenario\t\tRMS Attitude Error\tRMS Bias Error\tCoverage\n');
fprintf('--------\t\t------------------\t-------------\t--------\n');
for s = 1:length(scenarios)
    rms_att = rms(results{s}.attitude_error);
    rms_bias = rms(results{s}.bias_error);
    avg_cov = mean(results{s}.coverage);
    fprintf('%-15s\t%.6f\t\t%.2e\t\t%.1f%%\n', ...
        scenarios{s}.name, rms_att, rms_bias, avg_cov);
end

%% Main EKF function
function result = run_attitude_ekf(scenario, q_true, beta_true, omega_true, dt, N, sigma_star)
    % Initialize arrays
    time = (0:N-1)*dt;
    
    % True state propagation
    q_true_hist = zeros(4, N);
    beta_true_hist = zeros(3, N);
    q_true_hist(:,1) = q_true;
    beta_true_hist(:,1) = beta_true;
    
    % EKF initialization
    x_est = [q_true; zeros(3,1)]; % Perfect initial knowledge
    P = diag([0.1^2 * (pi/180)^2 * ones(1,3), (0.2 * pi/180/3600)^2 * ones(1,3)]);
    
    % Storage
    attitude_error = zeros(1, N);
    bias_error = zeros(1, N);
    coverage = zeros(1, N);
    
    % Reference vectors for star camera
    r_ref = eye(3);
    
    for k = 1:N
        % True state for comparison
        if k > 1
            [q_true_hist(:,k), beta_true_hist(:,k)] = propagate_true(...
                q_true_hist(:,k-1), beta_true_hist(:,k-1), omega_true, dt);
        end
        
        % Generate measurements
        gyro_meas = omega_true + beta_true_hist(:,k) + scenario.sigma_v * randn(3,1);
        star_meas = generate_star_measurement(q_true_hist(:,k), r_ref, sigma_star);
        
        % EKF Prediction
        [x_pred, P_pred] = ekf_predict(x_est, P, gyro_meas, scenario.sigma_u, scenario.sigma_v, dt);
        
        % EKF Update
        [x_est, P] = ekf_update(x_pred, P_pred, star_meas, r_ref, sigma_star);
        
        % Store results
        attitude_error(k) = quaternion_error(q_true_hist(:,k), x_est(1:4));
        bias_error(k) = norm(beta_true_hist(:,k) - x_est(5:7));
        
        % Check if error is within 3-sigma bounds
        sigma_bound = 3*sqrt(P(1,1)); % Simplified - use first attitude variance
        coverage(k) = (attitude_error(k) < sigma_bound) * 100;
    end
    
    result = struct();
    result.time = time;
    result.attitude_error = attitude_error;
    result.bias_error = bias_error;
    result.coverage = coverage;
end

%% Helper functions
function [q_new, beta_new] = propagate_true(q, beta, omega, dt)
    % Simple quaternion propagation
    omega_norm = norm(omega);
    if omega_norm > 1e-12
        axis = omega/omega_norm;
        angle = omega_norm * dt;
        dq = [axis * sin(angle/2); cos(angle/2)];
        q_new = quaternion_multiply(q, dq);
    else
        q_new = q;
    end
    q_new = q_new / norm(q_new);
    beta_new = beta; % Constant bias for simulation
end

function star_meas = generate_star_measurement(q_true, r_ref, sigma_star)
    % Star camera measurement (body vectors in reference frame)
    R = quaternion_to_rotation(q_true);
    star_meas = [];
    for i = 1:3
        vec = R(:,i) + sigma_star * randn(3,1);
        star_meas = [star_meas; vec/norm(vec)];
    end
end

function [x_pred, P_pred] = ekf_predict(x, P, gyro_meas, sigma_u, sigma_v, dt)
    % Simplified EKF prediction
    q = x(1:4); beta = x(5:7);
    
    % State prediction (simple rotation)
    omega_est = gyro_meas - beta;
    omega_norm = norm(omega_est);
    if omega_norm > 1e-12
        axis = omega_est/omega_norm;
        angle = omega_norm * dt;
        dq = [axis * sin(angle/2); cos(angle/2)];
        q_pred = quaternion_multiply(q, dq);
    else
        q_pred = q;
    end
    q_pred = q_pred / norm(q_pred);
    
    x_pred = [q_pred; beta];
    
    % Simplified covariance prediction
    F = eye(6); % Approximate state transition
    Q = diag([sigma_v^2 * dt * ones(1,3), sigma_u^2 * dt * ones(1,3)]);
    P_pred = F * P * F' + Q;
end

function [x_update, P_update] = ekf_update(x_pred, P_pred, star_meas, r_ref, sigma_star)
    % Simplified EKF update
    q_pred = x_pred(1:4);
    R_pred = quaternion_to_rotation(q_pred);
    
    % Measurement prediction
    h_pred = [];
    for i = 1:3
        h_pred = [h_pred; R_pred * r_ref(:,i)];
    end
    
    % Simplified measurement matrix (only use first vector for simplicity)
    H = zeros(3, 6);
    H(1:3, 1:3) = eye(3); % Only update attitude states
    
    % Kalman gain - FIXED: Use only first 3 measurements
    R_meas = sigma_star^2 * eye(3);
    K = P_pred * H' / (H * P_pred * H' + R_meas);
    
    % State update
    innovation = star_meas(1:3) - h_pred(1:3); % Use only first vector
    dx = K * innovation;
    
    x_update = x_pred;
    x_update(1:4) = quaternion_update(q_pred, dx(1:3));
    x_update(1:4) = x_update(1:4) / norm(x_update(1:4));
    
    % Covariance update
    P_update = (eye(6) - K * H) * P_pred;
end

function q_out = quaternion_multiply(q1, q2)
    % Quaternion multiplication
    q_out = [q1(4)*q2(1:3) + q2(4)*q1(1:3) + cross(q1(1:3), q2(1:3));
             q1(4)*q2(4) - dot(q1(1:3), q2(1:3))];
end

function R = quaternion_to_rotation(q)
    % Convert quaternion to rotation matrix
    q = q/norm(q);
    q0 = q(4); q1 = q(1); q2 = q(2); q3 = q(3);
    
    R = [1-2*(q2^2+q3^2), 2*(q1*q2-q0*q3), 2*(q1*q3+q0*q2);
         2*(q1*q2+q0*q3), 1-2*(q1^2+q3^2), 2*(q2*q3-q0*q1);
         2*(q1*q3-q0*q2), 2*(q2*q3+q0*q1), 1-2*(q1^2+q2^2)];
end

function q_new = quaternion_update(q, delta_theta)
    % Update quaternion with small rotation vector
    delta_q = [0.5 * delta_theta; 1];
    delta_q = delta_q / norm(delta_q);
    q_new = quaternion_multiply(q, delta_q);
end

function error = quaternion_error(q1, q2)
    % Calculate attitude error between two quaternions
    R1 = quaternion_to_rotation(q1);
    R2 = quaternion_to_rotation(q2);
    error_angle = acos((trace(R1'*R2) - 1)/2);
    error = abs(error_angle);
end