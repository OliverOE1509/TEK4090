% Exc 8

clear; close all; clc;
rng(42);

dt = 1; 
T = 60; 
N = T/dt;

omega_true = [0.001; 0.002; 0.0005]; 

r1 = [1;0;0]; r2 = [0;1;0]; r3 = [0;0;1];
r = [r1, r2, r3];

function S = skew_symmetric(v)
    S = [0, -v(3), v(2); v(3), 0, -v(1); -v(2), v(1), 0];
end

function q_out = quat_multiply(q1, q2)
    v1 = q1(1:3); s1 = q1(4);
    v2 = q2(1:3); s2 = q2(4);
    q_out = [s1*v2 + s2*v1 + cross(v1, v2); s1*s2 - dot(v1, v2)];
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

function euler = quat_to_euler(q)
    q0 = q(4); q1 = q(1); q2 = q(2); q3 = q(3);
    
    sinr_cosp = 2 * (q0 * q1 + q2 * q3);
    cosr_cosp = 1 - 2 * (q1 * q1 + q2 * q2);
    roll = atan2(sinr_cosp, cosr_cosp);
    sinp = 2 * (q0 * q2 - q3 * q1);
    if abs(sinp) >= 1
        pitch = sign(sinp) * pi/2;
    else
        pitch = asin(sinp);
    end
    
    siny_cosp = 2 * (q0 * q3 + q1 * q2);
    cosy_cosp = 1 - 2 * (q2 * q2 + q3 * q3);
    yaw = atan2(siny_cosp, cosy_cosp);
    
    euler = [roll; pitch; yaw];
end

function q_out = exact_attitude(q0, omega, t)
    omega_norm = norm(omega);
    if omega_norm > 1e-12
        axis = omega / omega_norm;
        angle = omega_norm * t;
        q_rot = [axis * sin(angle/2); cos(angle/2)];
        q_out = quat_multiply(q0, q_rot);
    else
        q_out = q0;
    end
    q_out = q_out / norm(q_out);
end

function [q_est, beta_est, P] = attitude_ekf_corrected(r, q0, beta0, P0, gyro_meas, star_meas, dt, sigma_v, sigma_u, sigma_star)
    r1 = r(:,1); r2 = r(:,2); r3 = r(:,3);
    N = size(gyro_meas, 2);
    q_est = zeros(4, N);
    beta_est = zeros(3, N);
    P = zeros(6, 6, N);
    
    q_est(:,1) = q0;
    beta_est(:,1) = beta0;
    P(:,:,1) = P0;
    
    for k = 1:N-1
        omega_est = gyro_meas(:,k) - beta_est(:,k);
        omega_norm = norm(omega_est);
        if omega_norm > 1e-12
            psi = sin(0.5 * omega_norm * dt) * omega_est / omega_norm;
            cos_term = cos(0.5 * omega_norm * dt);
            
            Omega = [cos_term * eye(3) - skew_symmetric(psi),  psi;
                     -psi',                                 cos_term];
            
            q_est(:,k+1) = Omega * q_est(:,k);
        else
            q_est(:,k+1) = q_est(:,k);
        end
        q_est(:,k+1) = q_est(:,k+1) / norm(q_est(:,k+1));
        if q_est(4,k+1) < 0
            q_est(:,k+1) = -q_est(:,k+1); 
        end
        beta_est(:,k+1) = beta_est(:,k);
        
        omega_norm = norm(omega_est);
        if omega_norm > 1e-12
            omega_skew = skew_symmetric(omega_est);
            Phi11 = eye(3) - omega_skew * sin(omega_norm * dt)/omega_norm + ...
                    omega_skew^2 * (1 - cos(omega_norm * dt))/(omega_norm^2);
            Phi12 = omega_skew * (1 - cos(omega_norm * dt))/(omega_norm^2) - ...
                    eye(3) * dt - omega_skew^2 * (omega_norm * dt - sin(omega_norm * dt))/(omega_norm^3);
        else
            Phi11 = eye(3);
            Phi12 = -eye(3) * dt;
        end
        Phi21 = zeros(3,3);
        Phi22 = eye(3);
        
        Phi = [Phi11, Phi12; Phi21, Phi22];
        
  
        Q11 = (sigma_v^2 * dt) * eye(3);
        Q12 = zeros(3,3);
        Q22 = (sigma_u^2 * dt) * eye(3);
        Qd = [Q11, Q12; Q12', Q22];
        
        P_pred = Phi * P(:,:,k) * Phi' + Qd;

        if ~isempty(star_meas) && k < size(star_meas, 2)
            q_pred = q_est(:,k+1);
            R_pred = quat_to_dcm(q_pred);
            
            H = zeros(9, 6);
            H(1:3, 1:3) = skew_symmetric(R_pred * r1);
            H(4:6, 1:3) = skew_symmetric(R_pred * r2);
            H(7:9, 1:3) = skew_symmetric(R_pred * r3);
            
            R_meas = (10 * sigma_star)^2 * eye(9);

            S = H * P_pred * H' + R_meas;
            reg_factor = max(1e-6, 1e-4 * dt) * eye(size(S));
            S = S + reg_factor;
            K = P_pred * H' / (S + 1e-8 * eye(size(S))); 
            
            b1_pred = R_pred * r1;
            b2_pred = R_pred * r2; 
            b3_pred = R_pred * r3;
            y_pred = [b1_pred; b2_pred; b3_pred];
            
            residual = star_meas(:,k+1) - y_pred;
            
            if norm(residual) > 10
                fprintf('k=%d: Large residual = %.3f, limiting...\n', k, norm(residual));
                residual = residual / norm(residual) * 0.1; 
            end
            delta_x = K * residual;
            
            max_att_update = 0.01 * sqrt(dt); 
            max_bias_update = 1e-8 * dt; 
            
            delta_x(1:3) = sign(delta_x(1:3)) .* min(abs(delta_x(1:3)), max_att_update);
            delta_x(4:6) = sign(delta_x(4:6)) .* min(abs(delta_x(4:6)), max_bias_update);
            
            delta_alpha = delta_x(1:3);
            delta_beta = delta_x(4:6);
            
            P_upd = (eye(6) - K * H) * P_pred;
            P_upd = 0.5 * (P_upd + P_upd');
            P_upd = P_upd + 1e-12 * eye(6);
            P(:,:,k+1) = 0.5 * (P_upd + P_upd');
            
            delta_q = [0.5 * delta_alpha; 1];
            delta_q = delta_q / norm(delta_q);
            q_est(:,k+1) = quat_multiply(delta_q, q_pred);
            q_est(:,k+1) = q_est(:,k+1) / norm(q_est(:,k+1));
            if q_est(4,k+1) < 0
                q_est(:,k+1) = -q_est(:,k+1);
            end

            
            beta_est(:,k+1) = beta_est(:,k+1) + delta_beta;
            
        else
            P(:,:,k+1) = P_pred;
        end
    end
end

scenarios = {
    struct('sigma_u', 1e-5, 'sigma_v', 1e-4, 'name', 'Baseline'),
    struct('sigma_u', 1e-5, 'sigma_v', 2e-4, 'name', 'High_sigma_v'), 
    struct('sigma_u', 2e-5, 'sigma_v', 1e-4, 'name', 'High_sigma_u')
};

sigma_star = 0.1 * pi/180;

results = cell(1,3);

for scenario_idx = 1:3 
    scenario = scenarios{scenario_idx};
    fprintf('\n=== Running scenario: %s ===\n', scenario.name);
    
    q_true_history = zeros(4, N);
    beta_true_history = zeros(3, N);
    euler_true_history = zeros(3, N);
    
    q_true_history(:,1) = [0; 0; 0; 1];
    true_beta_rad = [0.1; 0.1; 0.1] * pi/180/3600; 
    
    for k = 1:N
        if k == 1
            q_true_history(:,k) = [0; 0; 0; 1];
        else
            q_true_history(:,k) = exact_attitude(q_true_history(:,k-1), omega_true, dt);
        end
        beta_true_history(:,k) = true_beta_rad;
        euler_true_history(:,k) = quat_to_euler(q_true_history(:,k));
    end
    
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
    
    q0_est = q_true_history(:,1); 
    beta0_est = [0; 0; 0];

    P0a = (0.1 * pi/180)^2 * dt; 
    P0b = (0.1 * pi/180/3600)^2 * dt;   
    P0 = diag([P0a, P0a, P0a, P0b, P0b, P0b]);
    
    [q_est, beta_est, P] = attitude_ekf_corrected(r, q0_est, beta0_est, P0, gyro_meas, star_meas, dt, ...
                                       scenario.sigma_v, scenario.sigma_u, sigma_star);
    
    attitude_errors_deg = zeros(1, N);
    attitude_error_norms = zeros(1, N);
    euler_errors_deg = zeros(3, N);
    bias_errors_deg_hr = zeros(3, N);
    omega_errors_rad = zeros(3, N);
    omega_error_norms = zeros(1, N);
    
    for k = 1:N
        q_err = quat_multiply(quat_conjugate(q_est(:,k)), q_true_history(:,k));
        angle_err = 2 * acos(abs(q_err(4))) * 180/pi;
        attitude_errors_deg(k) = min(angle_err, 180-angle_err); 
        
        attitude_error_norms(k) = norm(q_err(1:3)) * 180/pi;
        
        euler_est = quat_to_euler(q_est(:,k));
        euler_errors_deg(:,k) = (euler_est - euler_true_history(:,k)) * 180/pi;
        
        bias_error_rad = beta_est(:,k) - beta_true_history(:,k);
        bias_errors_deg_hr(:,k) = bias_error_rad * (180/pi * 3600);
        omega_est = gyro_meas(:,k) - beta_est(:,k);
        omega_errors_rad(:,k) = omega_est - omega_true;
        omega_error_norms(k) = norm(omega_errors_rad(:,k)) * 180/pi; 
    end
    
    attitude_variances = zeros(3, N);
    bias_variances_deg_hr = zeros(3, N);
    
    for k = 1:N
        attitude_variances(:,k) = diag(P(1:3,1:3,k)) * (180/pi)^2;
        bias_variances_deg_hr(:,k) = diag(P(4:6,4:6,k)) * (180/pi * 3600)^2;
    end
    
    results{scenario_idx} = struct();
    results{scenario_idx}.name = scenario.name;
    results{scenario_idx}.attitude_errors_deg = attitude_errors_deg;
    results{scenario_idx}.attitude_error_norms = attitude_error_norms;
    results{scenario_idx}.euler_errors_deg = euler_errors_deg;
    results{scenario_idx}.bias_est_deg_hr = beta_est * (180/pi * 3600);
    results{scenario_idx}.bias_errors_deg_hr = bias_errors_deg_hr;
    results{scenario_idx}.attitude_variances = attitude_variances;
    results{scenario_idx}.bias_variances_deg_hr = bias_variances_deg_hr;
    results{scenario_idx}.P = P;
    results{scenario_idx}.omega_error_norms = omega_error_norms;
    results{scenario_idx}.omega_errors_rad = omega_errors_rad;
    
    fprintf('Final attitude error: %.4f deg\n', attitude_errors_deg(end));
    fprintf('Final attitude error norm: %.4f deg\n', attitude_error_norms(end));
    fprintf('Final bias errors: [%.3f, %.3f, %.3f] deg/hr\n', bias_errors_deg_hr(:,end));
end

time_minutes = (0:N-1)*dt/60;

colors = {'b', 'r', 'g'};
scenario_names = {scenarios{1}.name, scenarios{2}.name, scenarios{3}.name};

figure('Position', [100, 100, 1200, 800]);

subplot(2,2,1);
hold on;
for idx = 1:3
    plot(time_minutes, results{idx}.omega_error_norms, colors{idx}, 'LineWidth', 2, 'DisplayName', scenario_names{idx});
end
ylabel('Angular Velocity Error Norm (deg/s)');
xlabel('Time (min)');
title('Norms: estimated vs true angular velocities');
legend('show');
grid on;

subplot(2,2,2);
hold on;
for idx = 1:3
    bias_errors = results{idx}.bias_errors_deg_hr;
    bias_error_norms = sqrt(sum(bias_errors.^2, 1));
    plot(time_minutes, bias_error_norms, colors{idx}, 'LineWidth', 2, 'DisplayName', scenario_names{idx});
end
ylabel('Bias error borm (deg/hr)');
xlabel('Time (min)');
title('Norms: actual vs estimated gyro bias');
legend('show');
grid on;

subplot(2,2,3);
hold on;
for idx = 1:3
    plot(time_minutes, results{idx}.attitude_errors_deg, colors{idx}, 'LineWidth', 2, 'DisplayName', scenario_names{idx});
end
ylabel('Attitude Error (deg)');
xlabel('Time (min)');
title('Attitude norm comparison');
legend('show');
grid on;

subplot(2,2,4);
final_att_errors = [results{1}.attitude_errors_deg(end), results{2}.attitude_errors_deg(end), results{3}.attitude_errors_deg(end)];
final_bias_errors = [mean(abs(results{1}.bias_errors_deg_hr(:,end))), ...
                     mean(abs(results{2}.bias_errors_deg_hr(:,end))), ...
                     mean(abs(results{3}.bias_errors_deg_hr(:,end)))];
bar_data = [final_att_errors; final_bias_errors]';
bar(1:3, bar_data);
ylabel('Final Error');
xlabel('Scenario');
title('Final performance comparison');
legend('Attitude Error (deg)', 'Bias Error (deg/hr)', 'Location', 'best');
set(gca, 'XTickLabel', scenario_names);
grid on;