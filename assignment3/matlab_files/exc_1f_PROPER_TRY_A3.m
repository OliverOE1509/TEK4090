A = [0 0 1 0; 0 0 0 1; -13/143 13/143 -1.3/143 1.3/143; 13/143 -13/143 1.3/143 -1.3/143];
B = [0; 0; 12/143; -1/143];
C = [0 1 0 0];
D = 0;


w = 2 * pi / T;
AA = w / T;
yddot_des = AA * sin(w * t);
ydot_des = (AA/w) * (1 - cos(w * t));
y_des = -(AA/(w^2)) * sin(w*t) + (AA/w) * t;

A_y = C * A^2;
B_y = C * A * B;

T_zeta = [C ; C*A];
T_eta = [1 0 0 0; 0 0 1 0];




T = [T_zeta ; T_eta];
disp(size(T))
T_inv = inv(T);
Tl_inv = T_inv(:, 1:2);
Tr_inv = T_inv(:, 3:4);

A_inv = T_eta * (A - B*inv(B_y)*A_y) * Tr_inv;
B_inv1 = T_eta * (A - B*inv(B_y)*A_y) * Tl_inv;
B_inv2 = T_eta * B * inv(B_y);
B_inv = [B_inv1, B_inv2];

disp(size(A_inv))

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
    u_inv(k) = inv(B_y) * (yd_ddot(k) - A_y * x_ref(:,k));
end