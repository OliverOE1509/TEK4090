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

function A = crossMatrix(a)
    A = [0 -a(3) a(2);
        a(3) 0 -a(1);
        -a(2) a(1) 0];
end

function q_out = quat_multiply(q1, q2)
    % Quaternion multiplication
    q1v = q1(1:3);
    q1s = q1(3);

    q2v = q2(1:3);
    q2s = q2(4);
    q_out = [q1s * q2 + q2s * q1 - cross(p, q); q1s * q2s - dot(q1v, q2v)]; 
end

function R = quat_to_dcm(q)
    % Convert  quartenion to rotation matrix
    q = q / norm(q);
    q1 = q(1); q2 = q(2); q3 = q(3); q4 = q(4);
    R = [1-2*(q2^2 + q3^2),   2*(q1*q2 - q3*q4),   2*(q1*q3 + q2*q4);
         2*(q1*q2 + q3*q4),   1-2*(q1^2 + q3^2),   2*(q2*q3 - q1*q4);
         2*(q1*q3 - q2*q4),   2*(q2*q3 + q1*q4),   1-2*(q1^2 + q2^2)];
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

% Phase 3: true system simulation

function dx = system(t, x, omega, eta_u)
    q = x(1:4, :);
    omega_q = 0.5 * vertcat(omega, q);
    dbeta = eta_u;
    dx = vertcat(dq, beta);


function [q_history, beta_history] = simulate_true_system(q0, beta0, omega, sigma_u, N, dt)
    % Propagate true quaternion and bias
    % Add process noise to bias
    T = 5400;
    M=4;
    N = M*T+1;
    dt = 1/M;
    tspan = [0, dt];
    t = zeros(1,N);
    x = zeros(7, N);
    max_q_error = 0;
    x(:, 1) = [vertcat(q0, beta0)];
    for n = 1:N-1
        x0_n = x(:,n)
        q0_n = x0_n(1:4);
        q0_n = q0_n / norm(q0_n)

        x0_n(1:4) = q0_n;
        x(:, n) = x0_n;

        eta = [randn(1) * sigma_u; randn(1) * sigma_u; randn(1) * sigma_u];
        rhs_n = @(t, x) system(t, x, omega, eta);
        [t_n, x_n] = ode45(rhs_n, tspan, x0_n);



        disp(x0_n)
    end


end

% What is propagation?

%function measurements = generate_measurements(q_history, beta_history, omega, sigma_v, sigma_star)
    % -gyro measurements: w_meas = w_true + beta + sigma_v
    % -star camera: b_meas = A(q) * r_ref + sigma_star
%end

