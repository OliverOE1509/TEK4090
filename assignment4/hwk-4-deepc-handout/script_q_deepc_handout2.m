clear all; close all; clc;

addpath(genpath('./'))

% define system matrices and parameters

dt = 0.1;                                   % sampling time

A = [1     dt     0      0;
     0     1      0.05   0;
     0     -0.05  1      dt;
     0     0      0      1];
B = [0   0;
     dt  0;
     0   0;
     0   dt];
C = [1 0 0 0;
     0 0 1 0];
D = zeros(2,2);

[nz, nu] = size(B);                         % system size
[ny,  ~] = size(C);                         % output size

%% Offline data generation (for DeePC)
T_data = 400;                               % length of data
rng(0);
u_data = 0.5*randn(nu, T_data);             % persistently exciting input
z = zeros(nz,1);
y_data = zeros(ny, T_data);

for k = 1:T_data
    y_data(:,k) = C*z + D*u_data(:,k) + 0.001*randn([ny,1]);
    z = A*z + B*u_data(:,k);
end


%% DeePC parameters
T_ini = 10;                                 % length of "initial" window
n_pred= 20;                                 % prediction horizon
L     = T_ini + n_pred;

% Build Hankel matrices for U and Y
U_h = hankel_lti(u_data, L);                % size: (nu*L) x T
Y_h = hankel_lti(y_data, L);                % size: (ny*L) x T

T_cols = size(U_h,2);                       % number of columns = dim of g

% Partition into past (ini) and future (f)
U_p = U_h(1:nu*T_ini, :);                   % past inputs
U_f = U_h(nu*T_ini+1:end, :);               % future inputs
Y_p = Y_h(1:ny*T_ini, :);                   % past outputs
Y_f = Y_h(ny*T_ini+1:end, :);               % future outputs

% Cost and constraints
% Stage cost: sum_k (y_k - r_k)' Qy (y_k - r_k) + u_k' Ru u_k
Qy = diag([10, 10]);                % weight on outputs (2x2)
Ru = diag([1,  1]);                 % weight on inputs (2x2)

lambda_g = 1e-1;                    % regularization for g
lambda_s = 1e2;                     % penalty on slack s (y_ini constraint)

Qbar = kron(eye(n_pred), Qy);       % (ny*n_pred x ny*n_pred)
Rbar = kron(eye(n_pred), Ru);       % (nu*n_pred x nu*n_pred)

% Input constraints (same bounds for both inputs here)
u_min = -1.0;
u_max =  1.0;

% Output constraints (box constraints on each output component)
y_min = -0.5 * ones(ny,1);          % keep outputs within [-0.5, 1.5]
y_max =  1.5 * ones(ny,1);

% Closed-loop simulation with DeePC
n_sim = 80;                                 % simulation length
t_space = linspace(0,dt*n_sim,n_sim+1);     % simulation time space
x = [0; 0; 0; 0];                           % initial state (true system)
u_cl = zeros(nu, n_sim);
y_cl = zeros(ny, n_sim);

% reference output: step for both outputs
y_ref_step = [1; 1];                            % desired step
yd = 10*cos(2*pi*t_space/50);                   % desired trajectory
y_ref_vec  = repmat(y_ref_step, n_pred, 1);     % constant trajectory
%y_ref_vec = zeros(ny*n_pred, 1);
%for k = 1:n_pred
%    % Scale cosine to fit within [-0.5, 1.5] bounds
%    scaled_cos = 0.5 + cos(2*pi*(k)/50);  % Shifts from [-1,1] to [0,2], then scale
%   y_ref_vec((k-1)*ny+1:k*ny) = scaled_cos * [0.75; 0.75];  % Further scale to [0,1.5]
%end

% warm-up history (T_ini samples) with zero input
u_hist = zeros(nu, T_ini);
y_hist = zeros(ny, T_ini);

% quadprog options
opts = optimoptions('quadprog','Display','off');

% Dimensions
n_g = T_cols;               % dim of g
n_s = ny * T_ini;           % dim of slack s
n_z = n_g + n_s;            % total decision variable dim

% Pre-build constant parts of H and A_ineq for speed
H_gg = 2*( Y_f' * Qbar * Y_f + ...
           U_f' * Rbar * U_f + ...
           lambda_g * eye(n_g) );
H_ss = 2*lambda_s * eye(n_s);
H = blkdiag(H_gg, H_ss);

f_base = -2*( Y_f' * Qbar * y_ref_vec );  % part that does not depend on s
fs = zeros(n_s,1);
f = [f_base;
     fs];

% Input bounds
A_u = [ U_f  zeros(nu*n_pred, n_s);
       -U_f  zeros(nu*n_pred, n_s)];
b_u = [ u_max*ones(nu*n_pred,1);
       -u_min*ones(nu*n_pred,1)];

% Output bounds
A_y = [ Y_f  zeros(ny*n_pred, n_s);
       -Y_f  zeros(ny*n_pred, n_s)];
b_y = [ kron(ones(n_pred,1), y_max);
       -kron(ones(n_pred,1), y_min)];

% Slack positivity: s >= 0  <=>  -s <= 0
A_s = [zeros(n_s, n_g), -eye(n_s)];
b_s = zeros(n_s,1);

A_ineq = [A_u;
          A_y;
          A_s];
b_ineq = [b_u;
          b_y;
          b_s];

y_ref_full = zeros(ny, n_sim);
for k = 1:n_sim
    ct = t_space(k);
    cosVal = cos(2*pi*ct/50);
    scaled_value = 0.5 + 0.5*cosVal;
    y_ref_full(:, k) = scaled_value * [1; 1];
end

for k = 1:n_sim
    y_ref_vec = zeros(ny*n_pred, 1);    
    for idx = 1:n_pred
        ref_idx = min(k + idx - 1, n_sim);        
        y_ref_vec((idx-1)*ny+1:idx*ny) = y_ref_full(:, ref_idx);
    end
    f_base = -2*( Y_f' * Qbar * y_ref_vec );
    f = [f_base; zeros(n_s,1)];

    % Build u_ini and y_ini from history (stacked columns)
    u_ini = reshape(u_hist, nu*T_ini, 1);   % (nu*T_ini x 1)
    y_ini = reshape(y_hist, ny*T_ini, 1);   % (ny*T_ini x 1)

    % DeePC equalities:
    % [U_p  0] [g]   = [u_ini]
    % [Y_p -I] [s]   = [y_ini]
    Aeq = [U_p, zeros(nu*T_ini, n_s);
           Y_p, -eye(n_s)];
    beq = [u_ini;
           y_ini];

    % Solve QP for z = [g; s]
    H = (H + H')/2;
    z_opt = quadprog(H, f, A_ineq, b_ineq, Aeq, beq, [], [], [], opts);

    g_opt = z_opt(1:n_g);
    % s_opt = z_opt(n_g+1:end);  % available if you want to inspect slacks

    % Optimal future input sequence
    u_future = U_f * g_opt;          % size nu*n_pred x 1
    u_k = u_future(1:nu);            % apply first input (2x1)

    % Apply to true system
    y_k = C*x + D*u_k;
    x   = A*x + B*u_k;

    % Store closed-loop signals
    u_cl(:,k) = u_k;
    y_cl(:,k) = y_k;

    % Update history buffers (FIFO of length T_ini)
    if T_ini > 1
        u_hist = [u_hist(:,2:end), u_k];
        y_hist = [y_hist(:,2:end), y_k];
    else
        u_hist = u_k;
        y_hist = y_k;
    end
end

%% Plot results
t = t_space(1:end-1);

% Create the scaled cosine reference for plotting
y_ref_plot = zeros(ny, n_sim);
for k = 1:n_sim
    current_time = t_space(k);
    cos_value = cos(2*pi*current_time/50);
    scaled_value = 0.5 + 0.5*cos_value;  % ∈ [0,1]
    y_ref_plot(:, k) = scaled_value * [1; 1];
end

figure;
subplot(3,1,1);
plot(t, y_cl(1,:), 'LineWidth',1.5); hold on;
plot(t, y_ref_plot(1,:), 'r--', 'LineWidth',1.5);  % Plot cosine reference!
yline(y_min(1),':'); yline(y_max(1),':');
grid on;
xlabel('k'); ylabel('y_1(k)');
legend('y_1','ref_1 (cosine)','y_{1,min}','y_{1,max}');
title('DeePC closed-loop outputs with constraints - Cosine tracking');

subplot(3,1,2);
plot(t, y_cl(2,:), 'LineWidth',1.5); hold on;
plot(t, y_ref_plot(2,:), 'r--', 'LineWidth',1.5);  % Plot cosine reference!
yline(y_min(2),':'); yline(y_max(2),':');
grid on;
xlabel('k'); ylabel('y_2(k)');
legend('y_2','ref_2 (cosine)','y_{2,min}','y_{2,max}');

subplot(3,1,3);
plot(t, u_cl(1,:), 'LineWidth',1.5); hold on;
plot(t, u_cl(2,:), 'LineWidth',1.5);
yline(u_min,':'); yline(u_max,':');
grid on;
xlabel('k'); ylabel('u_i(k)');
legend('u_1','u_2','u_{min}','u_{max}');
title('Control inputs');

sgtitle('Cosine tracking with noise ε = 0.000001, λ_s = 1e4, λ_g = 1e-4');