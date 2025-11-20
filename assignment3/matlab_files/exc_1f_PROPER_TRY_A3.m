A = [0 0 1 0; 0 0 0 1; -13/143 13/143 -1.3/143 1.3/143; 13/143 -13/143 1.3/143 -1.3/143];
B = [0; 0; 12/143; -1/143];
C = [0 1 0 0];
D = 0;

Tt = 10;
t = 0:0.01:Tt;
w = 2 * pi / Tt;
AA = w / Tt;
ydd_dot = AA * sin(w * t);
yd_dot = (AA/w) * (1 - cos(w * t));
y_d = -(AA/(w^2)) * sin(w*t) + (AA/w) * t;

A_y = C * A^2;
B_y = C * A * B;

T_zeta = [C ; C*A];
T_eta = [1 0 0 0; 0 0 1 0];




T = [T_zeta ; T_eta];
T_inv = inv(T);
Tl_inv = T_inv(:, 1:2);
Tr_inv = T_inv(:, 3:4);

A_inv = T_eta * (A - B*inv(B_y)*A_y) * Tr_inv;
B_inv1 = T_eta * (A - B*inv(B_y)*A_y) * Tl_inv;
B_inv2 = T_eta * B * inv(B_y);
B_inv = [B_inv1, B_inv2];


function eta_ref = solve_internal_dynamics(A_inv, B_inv, y_d, yd_dot, yd_ddot, t)
    % Split A_inv into stable and unstable parts
    [V, D] = eig(A_inv);
    eigenvalues = diag(D);
    
    % Separate stable (Re<0) and unstable (Re>0) modes
    stable_idx = real(eigenvalues) < 0;
    unstable_idx = real(eigenvalues) > 0;
    
    % Transformation to decouple the system
    T_split = V;
    T_split_inv = inv(T_split);
    
    % Decoupled system matrices
    A_decoupled = T_split_inv * A_inv * T_split;  % Should be (nearly) diagonal
    B_decoupled = T_split_inv * B_inv;
    
    % Extract stable and unstable subsystems
    A_s = A_decoupled(stable_idx, stable_idx);
    B_s = B_decoupled(stable_idx, :);
    A_u = A_decoupled(unstable_idx, unstable_idx);  
    B_u = B_decoupled(unstable_idx, :);
    
    % Form the desired output vector Y_d(t)
    Y_d = [y_d; yd_dot; yd_ddot]';  % Each row is [y_d, yd_dot, yd_ddot] at time t
    
    % Initialize
    n_stable = sum(stable_idx);
    n_unstable = sum(unstable_idx);
    eta_s = zeros(n_stable, length(t));
    eta_u = zeros(n_unstable, length(t));
    
    % SOLVE STABLE PART FORWARD IN TIME (causal)
    for k = 2:length(t)
        dt = t(k) - t(k-1);
        % Forward Euler for stable part
        eta_s(:,k) = eta_s(:,k-1) + dt * (A_s * eta_s(:,k-1) + B_s * Y_d(k-1,:)');
    end
    
    % SOLVE UNSTABLE PART BACKWARD IN TIME (noncausal)
    for k = length(t)-1:-1:1
        dt = t(k+1) - t(k);
        % Backward Euler for unstable part (solving backward in time)
        eta_u(:,k) = eta_u(:,k+1) - dt * (A_u * eta_u(:,k+1) + B_u * Y_d(k+1,:)');
    end
    
    % Combine stable and unstable parts
    eta_decoupled = zeros(size(A_inv,1), length(t));
    eta_decoupled(stable_idx, :) = eta_s;
    eta_decoupled(unstable_idx, :) = eta_u;
    
    % Transform back to original eta coordinates
    eta_ref = T_split * eta_decoupled;
end

eta_ref = solve_internal_dynamics(A_inv, B_inv, y_d, yd_dot, ydd_dot, t);

for k = 1:length(t)
    zeta_d = [y_d(k); yd_dot(k)];
    x_ref(:,k) = T_inv * [zeta_d; eta_ref(:,k)];
end

for k = 1:length(t)
    u_inv(k) = inv(B_y) * (ydd_dot(k) - A_y * x_ref(:,k));
end





figure('Position', [100, 100, 1200, 800]);

% Plot 1: Desired trajectory
% subplot(3, 2, 1);
% plot(t, y_d, 'b', 'LineWidth', 2); hold on;
% plot(t, yd_dot, 'r', 'LineWidth', 2);
% plot(t, ydd_dot, 'g', 'LineWidth', 2);
% title('Desired Output Trajectory');
% xlabel('Time (s)');
% ylabel('Amplitude');
% legend('Position y_d', 'Velocity \dot{y}_d', 'Acceleration \ddot{y}_d', 'Location', 'best');
% grid on;

tracking_error = y_d - x_ref(2,:);
fprintf('Integral of y error')
disp(trapz(tracking_error))
fprintf('integral of control')
disp(trapz(u_inv))


figure;
subplot(3, 1, 1);
plot(t, y_d, 'b--', 'LineWidth', 2); hold on;
plot(t, x_ref(2,:), 'r-', 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Position');
legend('Desired y_d', 'Reference x_{2,ref}', 'Location', 'best');
grid on;

subplot(3, 1, 2)
plot(t, u_inv);
xlabel('Time (s)')
ylabel('u_inv(t)')
legend('control u')
grid on;

subplot(3, 1, 3)
plot(t, tracking_error);
ylabel('tracking error')
grid on;
return
% Plot 2: Internal states eta
% subplot(3, 3, 2);
% plot(t, eta_ref(1,:), 'b', 'LineWidth', 2); hold on;
% plot(t, eta_ref(2,:), 'r', 'LineWidth', 2);
% title('Internal States \eta_{ref}');
% xlabel('Time (s)');
% ylabel('Amplitude');
% legend('\eta_1', '\eta_2', 'Location', 'best');
% grid on;

% Plot 3: Inverse input
subplot(3, 2, 2);
plot(t, u_inv, 'm', 'LineWidth', 2);
title('Inverse Input u_{inv}');
xlabel('Time (s)');
ylabel('Force');
grid on;

% Plot 4: Reference states
% subplot(3, 2, 2);
% plot(t, x_ref(1,:), 'b', 'LineWidth', 2); hold on;
% plot(t, x_ref(2,:), 'r', 'LineWidth', 2);
% title('Reference States: Positions');
% xlabel('Time (s)');
% ylabel('Position');
% legend('x_{1,ref}', 'x_{2,ref}', 'Location', 'best');
% grid on;

% Plot 5: Reference velocities
% subplot(3, 3, 4);
% plot(t, x_ref(3,:), 'b', 'LineWidth', 2); hold on;
% plot(t, x_ref(4,:), 'r', 'LineWidth', 2);
% title('Reference States: Velocities');
% xlabel('Time (s)');
% ylabel('Velocity');
% legend('\dot{x}_{1,ref}', '\dot{x}_{2,ref}', 'Location', 'best');
% grid on;

% Plot 7: Compare desired vs reference output
subplot(3, 2, 5);
plot(t, y_d, 'b--', 'LineWidth', 2); hold on;
plot(t, x_ref(2,:), 'r-', 'LineWidth', 1);
title('Output Comparison');
xlabel('Time (s)');
ylabel('Position');
legend('Desired y_d', 'Reference x_{2,ref}', 'Location', 'best');
grid on;

% % Plot 8: Control effort components
% subplot(3, 3, 8);
% % Break down u_inv components for analysis
% u_components = zeros(3, length(t));
% for k = 1:length(t)
%     zeta_d = [y_d(k); yd_dot(k)];
%     u_components(1,k) = inv(B_y) * ydd_dot(k);
%     u_components(2,k) = -inv(B_y) * A_y * Tl_inv * zeta_d;
%     u_components(3,k) = -inv(B_y) * A_y * Tr_inv * eta_ref(:,k);
% end
% plot(t, u_components(1,:), 'g', 'LineWidth', 1); hold on;
% plot(t, u_components(2,:), 'b', 'LineWidth', 1);
% plot(t, u_components(3,:), 'r', 'LineWidth', 1);
% plot(t, u_inv, 'k--', 'LineWidth', 2);
% title('Control Effort Components');
% xlabel('Time (s)');
% ylabel('Force');
% legend('Accel term', 'State term', 'Internal term', 'Total u_{inv}', 'Location', 'best');
% grid on;

% Plot 9: Noncausal behavior indicator
% subplot(3, 3, 9);
% [V, D] = eig(A_inv);
% eigenvalues = diag(D);
% stable_idx = real(eigenvalues) < 0;
% unstable_idx = real(eigenvalues) > 0;
% 
% % Show which modes are solved forward/backward
% if any(unstable_idx)
%     plot(t, ones(size(t)), 'r-', 'LineWidth', 3); hold on;
%     text(mean(t), 1.1, 'UNSTABLE: Solved BACKWARD', 'HorizontalAlignment', 'center', 'Color', 'r');
% end
% if any(stable_idx)
%     plot(t, zeros(size(t)), 'g-', 'LineWidth', 3);
%     text(mean(t), -0.1, 'STABLE: Solved FORWARD', 'HorizontalAlignment', 'center', 'Color', 'g');
% end
% ylim([-0.5, 1.5]);
% title('Solution Direction');
% xlabel('Time (s)');
% yticks([]);

sgtitle('Feedforward Control Analysis for Nonminimum-Phase System');