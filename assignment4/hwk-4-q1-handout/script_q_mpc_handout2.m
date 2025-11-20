clear all; close all; clc;

addpath(genpath('./fcns'));

%% System matrices and parameters
dt = 0.1;                                   % sampling time

A = [1, dt;
     0, 1];

B = [0; dt];

C = [1, 0];                                 % y_k = C z_k = position

[nz, nu] = size(B);                         % system size
[ny,  ~] = size(C);                         % output size

G = 0.025*eye(nz);                          % process noise
K = 0.05*eye(ny);                           % measurement noise

z0 = [4;0];                                 % initial condition

%% MPC parameters
n_pred = 10;                                % prediction horizon
n_sim  = 100*n_pred;                        % simulation horizon
t_space = linspace(0,dt*n_sim,n_sim+1);     % simulation time grid
yd = 10*cos(2*pi*t_space/50);               % desired trajectory (output ref)

R = 0.1*eye(nu);                            % penalty on control
Qy = eye(ny);                               % penalty on output tracking error

%% SOFT state constraints (only on position z1)
pos_min = -5;
pos_max =  5;

% IMPORTANT: remove hard bounds on z1 here (we will encode them with slacks)
lbz = [-inf; -inf];                         % no hard bounds on state
ubz = [ inf;  inf];

lbu = [-1];                                 % lower bound on control
ubu = [ 1];                                 % upper bound on control

%% Dimensions for decision variable
n_zdec = nz*n_pred;     % z_1,...,z_T
n_udec = nu*n_pred;     % u_0,...,u_{T-1}
n_su   = n_pred;        % σ^u_1,...,σ^u_T (upper slack on position)
n_sl   = n_pred;        % σ^l_1,...,σ^l_T (lower slack on position)
Nx     = n_zdec + n_udec + n_su + n_sl;

%% Equality constraints (dynamics): Aeq * x = beq
Omega = sparse(getOmega(A,n_pred,'lti'));      % (n_zdec x nz)
Psi   = sparse(getPsi(A,B,n_pred,'lti'));      % (n_zdec x n_udec)

Aeq_zu = [-eye(n_zdec), Psi];                  % enforces z = Omega*z0 + Psi*u
Aeq    = [Aeq_zu, zeros(n_zdec, n_su + n_sl)]; % slacks do not enter dynamics
beq    = Omega*z0;

%% Inequality constraints with slacks: Aineq * x <= bineq
% soften:  pos_min - σ^l_k <= z1_k <= pos_max + σ^u_k
% equivalently
%   z1_k - σ^u_k <= pos_max
%  -z1_k - σ^l_k <= -pos_min

Aineq = zeros(2*n_pred, Nx);
bineq = zeros(2*n_pred, 1);

for k = 1:n_pred
    row_up  = 2*k - 1;   % upper inequality row index
    row_low = 2*k;       % lower inequality row index

    % index of position component z1_k in stacked z
    idx_z1 = (k-1)*nz + 1;

    % indices of slack variables in full x
    idx_su = n_zdec + n_udec + k;           % σ^u_k
    idx_sl = n_zdec + n_udec + n_su + k;    % σ^l_k

    % upper: z1_k - σ^u_k <= pos_max
    Aineq(row_up, idx_z1) =  1;
    Aineq(row_up, idx_su) = -1;
    bineq(row_up)         = pos_max;

    % lower: -z1_k - σ^l_k <= -pos_min
    Aineq(row_low, idx_z1) = -1;
    Aineq(row_low, idx_sl) = -1;
    bineq(row_low)         = -pos_min;
end

%% Simple box bounds for (z,u,σ)
lb_z      = kron(ones(n_pred,1), lbz);       % (n_zdec x 1)
ub_z      = kron(ones(n_pred,1), ubz);       % (n_zdec x 1)
lb_u      = kron(ones(n_pred,1), lbu);       % (n_udec x 1)
ub_u      = kron(ones(n_pred,1), ubu);       % (n_udec x 1)
lb_slack  = zeros(n_su + n_sl, 1);           % σ >= 0
ub_slack  = inf(n_su + n_sl, 1);             % no upper bound

lb = [lb_z; lb_u; lb_slack];
ub = [ub_z; ub_u; ub_slack];

%% Matrices for output-tracking cost
% stacked mapping from z to y: y_stack = (I⊗C) * z_stack
Cy_block = kron(speye(n_pred), C);           % (ny*n_pred x n_zdec)
Qy_block = kron(speye(n_pred), Qy);          % (ny*n_pred x ny*n_pred)

% Components used every step
Ky = Cy_block' * Qy_block * Cy_block;        % for H_zz
Ly = Cy_block' * Qy_block;                   % for f_z = -2*Ly*yref

% control cost part (constant)
Hu = 2*kron(speye(n_pred), R);

% slack cost part (constant, large weight)
w_slack = 1000;                              % choose large
Hslack  = 2*w_slack*eye(n_su + n_sl);

%% Simulation arrays
z_array = zeros(nz, n_sim+1);
z_array(:,1) = z0;
y_array = zeros(ny, n_sim+1);
y_array(:,1) = C*z0 + sqrtm(K)*randn([ny,1]);
u_array = zeros(nu, n_sim);

% optimizer options
opts.Display = 'off';

%% MPC loop
for kk = 1:n_sim

    % build reference over prediction horizon
    % indices kk+1 ... kk+n_pred in yd (length n_sim+1)
    idx_end = min(kk + n_pred, length(yd));
    yref_segment = yd(kk+1:idx_end);
    % if near the end, pad with last value
    if numel(yref_segment) < n_pred
        yref_segment = [yref_segment, repmat(yd(end), 1, n_pred - numel(yref_segment))];
    end
    yref_stack = yref_segment(:); % (ny*n_pred x 1), ny=1 here
    
    % build quadratic cost for this step

    Q_reg = 0.01 * eye(nz);
    Qreg_block = kron(eye(n_pred), Q_reg);

    H_zz = 2*Ky + 2*Qreg_block;
    f_z  = -2*Ly*yref_stack;      % output tracking term

    H = blkdiag(H_zz, Hu, Hslack);
    f = [f_z;
         zeros(n_udec,1);
         zeros(n_su + n_sl,1)];

    % solve QP
    [x,~,exitflag] = quadprog(H,f,Aineq,bineq,Aeq,beq,lb,ub,[],opts);

    if exitflag <= 0 || isempty(x)
        fprintf('QP failed at iteration %d with exitflag %d\n', kk, exitflag);
        break
    end

    % unpack optimal solution
    z_opt = reshape(x(1:n_zdec), nz, n_pred);
    u_opt = reshape(x(n_zdec+1:n_zdec+n_udec), nu, n_pred);
    u_array(:,kk) = u_opt(:,1);

    % propagate and measure
    z_array(:,kk+1) = A*z_array(:,kk) + B*u_array(:,kk) + sqrtm(G)*randn([nz,1]);
    y_array(:,kk+1) = C*z_array(:,kk) + sqrtm(K)*randn([ny,1]);

    % update equality constraint rhs for next step
    beq = Omega*z_array(:,kk+1);
end

kk_end = kk;  % last successful step (or n_sim if finished)

%% Plots
figure;
subplot(3,1,1)
plot(t_space(1:kk_end), y_array(1,1:kk_end), 'b', ...
     t_space(1:kk_end), yd(1:kk_end),       'r--');
ylabel('y = position')
title('MPC with output tracking and slack variables')
legend('y','y_d')

subplot(3,1,2)
plot(t_space(1:kk_end), z_array(2,1:kk_end));
ylabel('z_2 (velocity)')

subplot(3,1,3)
plot(t_space(1:kk_end-1), u_array(1,1:kk_end-1));
ylabel('u (input)')
xlabel('time')
