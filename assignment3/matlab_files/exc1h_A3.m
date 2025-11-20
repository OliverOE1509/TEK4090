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

Q = diag([1, 10, 0.1, 0.1]);
R = 0.01;

K = lqr(A, B, Q, R);


% Heres the forward analysis solution of control u
x0 = [0; 0; 0; 0];  % Initial state
x = x0;
X_sim = []; U_fb = []; U_total = []; Y_sim = [];

for k=1:length(t)
    x_refd = x_ref(:, k);
    u_fb = -K * (x - x_refd);
    u_total = u_fb + u_inv(k);

    X_sim = [X_sim, x];
    U_fb = [U_fb, u_fb];
    U_total = [U_total, u_total];
    Y_sim = [Y_sim, x(2)];

    if k < length(t)
        dt = t(k+1) - t(k);
        xdot = A*x + B* u_total;
        x = x + xdot * dt;
    end
end

disp(size(U_fb))
disp(size(U_total))
disp(size(u_inv))

tracking_error = y_d - Y_sim;

fprintf('Integral of y error')
disp(trapz(tracking_error))
fprintf('integral of control')
disp(trapz(U_total))

figure('Position', [100, 100, 1200, 800]);

figure;
subplot(3, 1, 1);
plot(t, y_d, 'b--', 'LineWidth', 2); hold on;
plot(t, Y_sim, 'r-', 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Position');
legend('Desired y_d', 'Reference x_{2,ref}', 'Location', 'best');
grid on;

subplot(3, 1, 2)
plot(t, u_inv, 'g-', 'LineWidth', 2, 'DisplayName', 'Feedforward'); hold on;
plot(t, U_fb, 'r-', 'LineWidth', 1, 'DisplayName', 'Feedback');
plot(t, U_total, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Combined');xlabel('Time (s)')
ylabel('u_inv(t)')
legend('control u')
grid on;

subplot(3, 1, 3)

plot(t, tracking_error, 'r', 'LineWidth', 2);
ylabel('tracking error')
grid on;

