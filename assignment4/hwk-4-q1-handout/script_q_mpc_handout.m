clear all; close all; clc;

addpath(genpath('./fcns'));

% define system matrices and parameters

dt = 0.1;                                   % sampling time

A = [1, dt;
     0, 1];

B = [0; dt];

C = [1, 0];

[nz, nu] = size(B);                         % system size
[ny,  ~] = size(C);                         % output size

G = 0.025*eye(nz);                          % process noise
K = 0.05*eye(ny);                           % measurement noise

z0 = [4;0];                                 % initial condition

% define MPC parameters

n_pred = 10;                                % prediction horizon
n_sim  = 100*n_pred;                        % simulation horizon
t_space = linspace(0,dt*n_sim,n_sim+1);     % simulation time grid
yd = 10*cos(2*pi*t_space/50);               % desired trajectory

Q = eye(nz);                                % penalty on state variable
R = 0.1*eye(nu);                            % penalty on control variable

Qy = eye(ny);                               % penalty on output variable

n_zdec = nz * n_pred;
n_udec = nu * n_pred;
n_su = n_pred;
n_sl = n_pred;
Nx = n_zdec + n_udec + n_su + n_sl;

% define constraint matrices
lbz = [-inf;-inf];                            % lower bound on state
ubz = [inf; inf];                             % upper bound on state

lbu = [-1];                                 % lower bound on control
ubu = [1];                                  % upper bound on control

lb_z = kron(ones(n_pred, 1), lbz);
ub_z = kron(ones(n_pred, 1), ubz);
lb_u = kron(ones(n_pred, 1), lbu);
ub_u = kron(ones(n_pred, 1), ubu);

lb_slack = zeros(n_su + n_sl, 1);
ub_slack = inf(n_su + n_sl, 1); 

Omega = sparse(getOmega(A,n_pred,'lti'));
Psi   = sparse(getPsi(A,B,n_pred,'lti'));
Aeq_zu = [eye(n_zdec), -Psi];
Aeq = [Aeq_zu, zeros(n_zdec, n_su + n_sl)];
beq = Omega*z0;

Aineq = zeros(2*n_pred, Nx);
bineq = zeros(2*n_pred, 1);

pos_max = 5;
pos_min = -5;

for k = 1:n_pred
    row_up  = 2*k - 1;
    row_low = 2*k;

    % index of z1_k in stacked z
    idx_z1 = (k-1)*nz + 1;

    % indices of sigma^u_k and sigma^l_k in the full x
    idx_su = n_zdec + n_udec + k;
    idx_sl = n_zdec + n_udec + n_su + k;

    % z1_k - sigma^u_k <= pos_max
    Aineq(row_up, idx_z1) =  1;
    Aineq(row_up, idx_su) = -1;
    bineq(row_up)         =  pos_max;

    % -z1_k - sigma^l_k <= -pos_min
    Aineq(row_low, idx_z1) = -1;
    Aineq(row_low, idx_sl) = -1;
    bineq(row_low)         = -pos_min;   % = 5
end


% simulate MPC 

z_array = zeros(nz,n_sim);
z_array(:,1) = z0;
y_array = zeros(ny,n_sim);
y_array(:,1) = C*z0 + sqrtm(K)*randn([ny,1]);
u_array = zeros(nu,n_sim);

% optimizer options
opts.Display = 'off';

exitflags = zeros(1,n_sim);
kk_crash = NaN;

lb = [lb_z; lb_u; lb_slack];
ub = [ub_z; ub_u; ub_slack];

for kk = 1:n_sim
    fprintf('Iteration %i\n',kk)
    disp(z_array(:,kk))

    % solve MPC

    % 1) reference segment over prediction horizon
    idx_end = min(kk + n_pred, length(yd));
    yref_seg = yd(kk+1:idx_end);
    if numel(yref_seg) < n_pred
        yref_seg = [yref_seg, repmat(yd(end), 1, n_pred-numel(yref_seg))];
    end
    yref_stack = yref_seg(:);       % (ny*n_pred x 1), ny = 1 here
    
    % 2) matrices for output tracking: y_stack = (I⊗C) z_stack
    Cy_block = kron(speye(n_pred), C);        % (ny*n_pred x n_zdec)
    Qy_block = kron(speye(n_pred), Qy);       % (ny*n_pred x ny*n_pred)
    
    % tracking quadratic and linear parts in z
    Ky = Cy_block' * Qy_block * Cy_block;     % (n_zdec x n_zdec)
    Ly = Cy_block' * Qy_block;                % (n_zdec x ny*n_pred)
    
    % 3) original state cost Q and input cost R
    Q_block = kron(speye(n_pred), Q);         % from part (b)
    R_block = kron(speye(n_pred), R);
    
    % 4) slack cost
    w_slack = 1000;
    Hslack  = 2*w_slack*eye(n_su + n_sl);
    
    % 5) assemble H and f  (J = 1/2 x'Hx + f'x)
    H_zz = 2*(Q_block + Ky);                  % combine state and output terms
    f_z  = -2*Ly*yref_stack;                  % linear part from tracking
    
    Hu = 2*R_block;
    
    H = blkdiag(H_zz, Hu, Hslack);
    f = [f_z;
         zeros(n_udec,1);
         zeros(n_su + n_sl,1)];

    [x,~,exitflag] = quadprog(H,f,Aineq,bineq,Aeq,beq,lb,ub,[],opts);
    % detect infeasible QP
    if exitflag == -2 || isempty(x)
        kk_crash = kk;
        disp('QP infeasible. Stopping simulation.')
        disp('State that caused infeasibility:')
        disp(z_array(:,kk))
        break
    end

    exitflags(kk) = exitflag;

    % unpack optimal solution
    z_opt = reshape(x(1:n_zdec), nz, n_pred);
    u_opt = reshape(x(n_zdec+1:n_zdec+n_udec), nu, n_pred);
    u_array(:,kk) = u_opt(:,1);

    % propagate
    z_array(:,kk+1) = A*z_array(:,kk) + B*u_array(:,kk) + sqrtm(G)*randn([nz,1]);
    y_array(:,kk+1) = C*z_array(:,kk) + sqrtm(K)*randn([ny,1]);

    % update eq constraints
    beq = Omega*z_array(:,kk+1);
end


if isnan(kk_crash)
    kk_end = n_sim;
else
    kk_end = kk_crash;
end

figure;
subplot(3,1,1)
plot(t_space(1:kk_end+1), y_array(1,1:kk_end+1), 'b', ...
     t_space(1:kk_end+1), yd(1:kk_end+1), 'r--');
ylabel('y (position)')
legend('y','y_d')
title('MPC with tracking and slack')


subplot(3,1,2)
plot(t_space(1:kk_end), z_array(2,1:kk_end))
ylabel('z_2 (velocity)')

subplot(3,1,3)
plot(t_space(1:kk_end), u_array(1,1:kk_end))
ylabel('u (input)')
xlabel('time')