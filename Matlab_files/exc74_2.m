%% Q8
%Define system, parameters and functions

%Initial conditions
q_0 = [0;0;0;1]; %Initial attitude (quaternion)
beta_0 = [0.1;0.1;0.1]; %Initial bias

%Parameters
sigma_i = 0.01;   %Star camera noise
var_i = (sigma_i^2)*eye(3);
R = blkdiag(var_i, var_i, var_i);
omega = [0;0;pi/2700];   %Angular velocity
r1 = [1;0;0];    %Reference basis vector
r2 = [0;1;0];    %Reference basis vector
r3 = [0;0;1];    %Reference basis vector
r = [r1;r2;r3];

function A = crossMatrix(a)
A = [    0    -a(3)    a(2);
    a(3)    0    -a(1);
    -a(2)    a(1)    0 ];
end

function pq = psi_q(quat)
up = quat(4)*eye(3) - crossMatrix(quat(1:3));
down = (-1)*quat(1:3).';
pq = vertcat(up,down);
end

function pq = xi_q(quat)
up = quat(4)*eye(3) + crossMatrix(quat(1:3));
down = (-1)*quat(1:3).';
pq = vertcat(up,down);
end

function A = q_to_matrix(quat)
A = (xi_q(quat).')*(psi_q(quat));
end

function quat = q_times(q1,q2)
%{
Implements: quat = q1 \otimes q2
Args:
q1, q2 - Two quaternions. I.e. vertical, two dimensional vectors
of length 1.
Returns:
quat: The resulting quaternion
%}

quat = [psi_q(q1) q1]*q2;
end

function dx = rhs(t, x, omega, eta_u)
q = x(1:4,:);
omega_q = 0.5*vertcat(omega,0);
dq = q_times(omega_q, q);

dbeta = eta_u;

dx = vertcat(dq, dbeta);
end

%Scenario 1
%Gyro noise
s1_sigma_u = 1; %Bias drift
s1_sigma_v = 1; %Measurement noise
%Scenario 2
%Gyro noise
s2_sigma_u = 1; %Bias drift
s2_sigma_v = 2; %Measurement noise
%Scenario 3
%Gyro noise
s3_sigma_u = 2; %Bias drift
s3_sigma_v = 1; %Measurement noise

%Simulate true system
function [t, x, max_q_err] = simulate_sys(q_0, beta_0, omega, sigma_u)
T = 5400;
M = 4;
N = T*M+1;
dt = 1/M;
tspan = [0,dt];

t = zeros(1,N);
x = zeros(7,N);
max_q_err = 0;
x(:,1) = [vertcat(q_0, beta_0)];

for n = 1:N-1
    %Get state at the start of time interval
    x0_n = x(:,n);

    %Normalize quanterion
    q0_n = x0_n(1:4);
    q0_n = q0_n/norm(q0_n);
    x0_n(1:4) = q0_n;

    %Replace state
    x(:,n) = x0_n;

    %Draw eta
    eta1 = randn(1)*sigma_u;
    eta2 = randn(1)*sigma_u;
    eta3 = randn(1)*sigma_u;
    eta = [eta1;eta2;eta3];

    %Simulate interval
    rhs_n = @(t,x) rhs(t, x, omega, eta);
    [t_n, x_n] = ode45(rhs_n, tspan, x0_n);

    %Log max quanterion error
    q_norms = vecnorm(x_n(:,1:4), 2, 2); % Euclidean norm of each quaternion
    q_err = max(abs(1-q_norms));
    max_q_err = max([max_q_err, q_err]);

    %Log values at the end of the interval
    t(n+1) = t(n) + dt;
    x(:,n+1) = x_n(end,:).';

end
end

[s1_t, s1_x, s1_max_q_err] = simulate_sys(q_0, beta_0, omega, s1_sigma_u);
[s2_t, s2_x, s2_max_q_err] = simulate_sys(q_0, beta_0, omega, s2_sigma_u);
[s3_t, s3_x, s3_max_q_err] = simulate_sys(q_0, beta_0, omega, s3_sigma_u);

%Ensure they have the same true orientation
%I only want the biases to be different.
s2_x(1:4,:) = s1_x(1:4,:);
s3_x(1:4,:) = s1_x(1:4,:);

%Simulate measurement model
tol = 1e-9;
t_k = abs(s1_t - round(s1_t)) < tol;
s1_x_k = s1_x(:, t_k);
s2_x_k = s2_x(:, t_k);
s3_x_k = s3_x(:, t_k);

function bi_k = draw_body_vec(quat_k, ri, sigma_i)
v1 = randn(1)*sigma_i;
v2 = randn(1)*sigma_i;
v3 = randn(1)*sigma_i;
zi = [v1;v2;v3];

bi = q_to_matrix(quat_k)*ri + zi;
bi_k = bi/norm(bi);
end

function omega_k = draw_omega(omega, beta_k, sigma_v)
v1 = randn(1)*sigma_v;
v2 = randn(1)*sigma_v;
v3 = randn(1)*sigma_v;
eta_v = [v1;v2;v3];

omega_k = omega + beta_k + eta_v;
end

function s_y_k = simulate_meas(s_x_k, sigma_v, sigma_i, r1, r2, r3, omega)
J = length(s_x_k);
s_y_k = zeros(12, J);
for k = 1:J
    x_k = s_x_k(:,k);
    quat_k = x_k(1:4,:);
    beta_k = x_k(5:7,:);

    b1_k = draw_body_vec(quat_k, r1, sigma_i);
    b2_k = draw_body_vec(quat_k, r2, sigma_i);
    b3_k = draw_body_vec(quat_k, r3, sigma_i);

    omega_k = draw_omega(omega, beta_k, sigma_v);

    s_y_k(:,k) = vertcat(b1_k, b2_k, b3_k, omega_k);
end
end

s1_y_k = simulate_meas(s1_x_k, s1_sigma_v, sigma_i, r1, r2, r3, omega);
s2_y_k = simulate_meas(s2_x_k, s2_sigma_v, sigma_i, r1, r2, r3, omega);
s3_y_k = simulate_meas(s3_x_k, s3_sigma_v, sigma_i, r1, r2, r3, omega);

%I want them to have the same attitude measurements but different gyro
%measurements
s2_y_k(1:9,:) = s1_y_k(1:9,:);
s3_y_k(1:9,:) = s1_y_k(1:9,:);

%Init filter params
p0_a = 2; %Error covariance of Euler angles
p0_b = 2; %Error covariance of Gyro bias
p0 = diag([p0_a p0_a p0_a p0_b p0_b p0_b]);
e_q0 = [0; 0; 0.976327476659339; -0.216297615150996]; %Expected Initial attitude (quanternion)
e_beta0 = [1.020260509045736; -0.042714399068091; -0.233307253785185]; %Expected Initial gyro bias

%Implement Filter
function Hk = meas_jacobian(qk_m, r1, r2, r3)
    A = q_to_matrix(qk_m);
    Hk = [
    crossMatrix(A*r1) zeros(3,3);
    crossMatrix(A*r2) zeros(3,3);
    crossMatrix(A*r3) zeros(3,3)
    ];
end

function hk = meas_func(qk_m, r1, r2, r3)
    A = q_to_matrix(qk_m);
    hk = [A*r1; A*r2; A*r3];
end

function Kk = gain_k(pk_m, Hk, R)
    Kk = pk_m*(Hk.')*inv(Hk*pk_m*(Hk.') + R);
end

function [qk_p, betak_p, omegak_p, pk_p] = output_to_est(r, yk, qk_m, betak_m, pk_m, R)
    %Takes output at time k and updates.
    %Unpack and prepare variables
    bk_m = yk(1:9);
    omegak_m = yk(10:12);
    r1 = r(1:3);
    r2 = r(4:6);
    r3 = r(7:9);
    %Generate matrices.
    Hk = meas_jacobian(qk_m, r1, r2, r3);
    hk = meas_func(qk_m, r1, r2, r3);
    Kk = gain_k(pk_m, Hk, R);
    %Generate error state and unpack delta alpha and beta
    err_st = Kk*(bk_m-hk);
    delta_alpha = err_st(1:3);
    delta_beta = err_st(4:6);
    %Update covariance
    pk_p = (eye(6)-Kk*Hk)*pk_m;
    %Update quaternion
    qk_p = qk_m + 0.5*xi_q(qk_m)*delta_alpha;
    qk_p = qk_p /norm(qk_p);
    %Update beta
    betak_p = betak_m + delta_beta;
    %Update omega
    omegak_p = omegak_m - betak_p;
    %Return
end

function psik = psi_vector(omegak_p)
sclr = sin(0.5*norm(omegak_p))/norm(omegak_p);
psik = sclr*omegak_p;
end

function O = prop_q_matrix(omegak_p)
psik = psi_vector(omegak_p);
o_norm = norm(omegak_p);

O11 = cos(0.5*o_norm)*eye(3) - crossMatrix(psik);
O12 = psik;
O21 = (-1)*psik.';
O22 = cos(0.5*o_norm);

O = [O11 O12; O21 O22];
end

function gammak = prop_gamma_matrix()
gammak = [(-1)*eye(3) zeros(3,3); zeros(3,3) eye(3)];
end

function Qk = prop_Q_matrix(sigma_v, sigma_u)
    var_v = sigma_v^2;
    var_u = sigma_u^2;
    Q11 = ( var_v + (1/3)*var_u )*eye(3);
    Q12 = ( (1/2)*var_u )*eye(3);
    Q21 = ( (1/2)*var_u )*eye(3);
    Q22 = ( var_u )*eye(3);
    Qk = [Q11 Q12; Q21 Q22];
end

function thethak = prop_thetha_matrix(omegak_p)
    cmo = crossMatrix(omegak_p);
    cmo_2 = cmo^2;
    o_norm = norm(omegak_p);
    o_norm_2 = o_norm^2;

    thet11 = eye(3) - cmo*(sin(o_norm)/o_norm) + cmo_2*(-(1 - cos(o_norm))/(o_norm_2));
    thet12 = cmo*((1 - cos(o_norm))/(o_norm_2)) - eye(3) - cmo_2*(-(o_norm - sin(o_norm))/(o_norm^3));

    thet21 = zeros(3,3);
    thet22 = eye(3);
    thethak = [thet11 thet12; thet21 thet22];
end

function [qkl_m, betakl_m, pkl_m] = prop_est(qk_p, betak_p, omegak_p, pk_p, sigma_v, sigma_u)
    %Takes updated estimates and proposates
    %Proposate beta
    betakl_m = betak_p;
    %Proposate quaternion
    O = prop_q_matrix(omegak_p);
    qkl_m = O*qk_p;

    %Proposate covariance
    gammak = prop_gamma_matrix();
    Qk = prop_Q_matrix(sigma_v, sigma_u);
    thethak = prop_thetha_matrix(omegak_p);
    pkl_m = thethak*pk_p*(thethak.') + gammak*Qk*(gammak.');
end

%Simulate estimator
function [est_state, err_covar] = simulate_filter(r, R, sigma_u, sigma_v, y, p0, e_q0, e_beta0 )
    %Takes entire meas sequence + params and returns estimate sequence + errors
    N = length(y);
    est_state = zeros(10,N); %First quaternion, then bias, then angular velocity
    err_covar = zeros(6,6,N);

    %Est initial
    [qk_p, betak_p, omegak_p, pk_p] = output_to_est(r, y(:,1), e_q0, e_beta0, p0, R);
    est_state(1:4,1) = qk_p;
    est_state(5:7,1) = betak_p;
    est_state(8:10,1) = omegak_p;
    err_covar(:,:,1) = pk_p;
    %Est and prop open interval for k = 1:N-1
    for k = 1:N-1
        %Propogate
        qk_p = est_state(1:4,k);
        betak_p = est_state(5:7,k);
        omegak_p = est_state(8:10,k);
        pk_p = err_covar(:,:,k);
        [qkl_m, betakl_m, pkl_m] = prop_est(qk_p, betak_p, omegak_p, pk_p, sigma_v, sigma_u);

        %Update
        ykl = y(:,k+1);
        [qkl_p, betakl_p, omegakl_p, pkl_p] = output_to_est(r, ykl, qkl_m, betakl_m, pkl_m, R);

        est_state(1:4,k+1) = qkl_p;
        est_state(5:7,k+1) = betakl_p;
        est_state(8:10,k+1) = omegakl_p;
        err_covar(:,:,k+1) = pkl_p;
    end
    %Return
end

[s1_est_state, s1_err_covar] = simulate_filter(r, R, s1_sigma_u, s1_sigma_v, s1_y_k, p0, e_q0, e_beta0);
[s2_est_state, s2_err_covar] = simulate_filter(r, R, s2_sigma_u, s2_sigma_v, s2_y_k, p0, e_q0, e_beta0);
[s3_est_state, s3_err_covar] = simulate_filter(r, R, s3_sigma_u, s3_sigma_v, s3_y_k, p0, e_q0, e_beta0);

%Assess performance
function [ate, bie, ave] = calc_error(true_x, true_omega, est_x)
N = length(true_x);

ate = zeros(1,N);
bie = zeros(1,N);
ave = zeros(1,N);
for k = 1:N
    %Attitude error (Forbenius norm between attitude matrices)
    truA = q_to_matrix(true_x(1:4,k));
    estA = q_to_matrix(est_x(1:4,k));
    ate(k) = norm(truA-estA, 'fro');
    %Bias error (Euclidean norm between bias vectors)
    truBeta = true_x(5:7,k);
    estBeta = est_x(5:7,k);
    bie(k) = norm(truBeta-estBeta);
    %Angular velocity error (Euclidean norm between velocity vectors)
    estOmega = est_x(8:10, k);
    ave(k) = norm(true_omega-estOmega);
end
end

[s1_ate, s1_bie, s1_ave] = calc_error(s1_x_k, omega, s1_est_state);
[s2_ate, s2_bie, s2_ave] = calc_error(s2_x_k, omega, s2_est_state);
[s3_ate, s3_bie, s3_ave] = calc_error(s3_x_k, omega, s3_est_state);

function [avg_ea, avg_bia] = calc_avg_var(err_covar)
N = length(err_covar);

avg_ea = zeros(1,N);
avg_bia = zeros(1,N);
for k = 1:N
    elem = diag(err_covar(:,:,k));
    ea = elem(1:3);
    bia = elem(4:6);
    %Average variance of euler angles
    avg_ea(k) = mean(ea);
    %Average variance of Bias
    avg_bia(k) = mean(bia);
end
end

[s1_avg_ea, s1_avg_bia] = calc_avg_var(s1_err_covar);
[s2_avg_ea, s2_avg_bia] = calc_avg_var(s2_err_covar);
[s3_avg_ea, s3_avg_bia] = calc_avg_var(s3_err_covar);

%plot results

%1:
figure;
N = length(s1_x_k);
%Attitude error
subplot(3,1,1); plot(1:N, [s1_ate.',s2_ate.',s3_ate.']); title('Attitude error (Forbenius norm)'); legend('Scenario 1','Scenario 2','Scenario 3'); xlim([1, N]);
%Bias error
subplot(3,1,2); plot(1:N, [s1_bie.',s2_bie.',s3_bie.']); title('Gyro Bias error (Euclidean norm)'); legend('Scenario 1','Scenario 2','Scenario 3'); xlim([1, N]);
%Angular velocity error
subplot(3,1,3); plot(1:N, [s1_ave.',s2_ave.',s3_ave.']); title('Angular velocity error (Euclidean norm)'); legend('Scenario 1','Scenario 2','Scenario 3'); xlim([1, N]);

%2:
figure;
%Avg euler angle variance
subplot(2,1,1); plot(1:N, [s1_avg_ea.',s2_avg_ea.', s3_avg_ea.']); title('Average Variance of Euler Angle Residuals'); legend('Scenario 1','Scenario 2','Scenario 3'); xlim([1, N]);
%Avg bias variance
subplot(2,1,2); plot(1:N, [s1_avg_bia.',s2_avg_bia.', s3_avg_bia.']); title('Average Variance of Gyro Bias Residuals'); legend('Scenario 1','Scenario 2','Scenario 3'); xlim([1, N]);