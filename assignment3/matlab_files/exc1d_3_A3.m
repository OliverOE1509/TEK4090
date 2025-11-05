A_r = 1;
L_r = 1;
E_r = 1;
rao_r = 6;
M1 = 10;
M2 = 10;
alpha = 0.1;

M_lump = [10 0; 0 10];
M_rod = (A_r * L_r * rao_r) / 6 * [2 1; 1 2];
K_rod = (A_r * E_r) / L_r * [1 -1; -1 1];
C_rod = alpha * K_rod;

A = [0 0 1 0;
    0 0 0 1;
    -13/143 13/143 -1.3/143 1.3/143;
    13/143 -13/143 1.3/143 -1.3/143];
B = [0;0;12/143;-1/143];
C = [0 1 0 0];
D = 0;




O = obsv(A, C);

sys = ss(A, B, C, D);
poles = eig(A);
zeros = tzero(sys);

T = 10;
w = 2 * pi / T;
A_a = w / T;
t_ref = 0:0.01:T;



Q = diag([1, 1000, 10, 100]);

A_aug = [A, ones(4,1) * 0; -C, 0];
B_aug = [B; 0];
Q_aug = diag([1, 1000, 10, 100, 1000]);
R = 0.1;
K = lqr(A_aug, B_aug, Q_aug, R);

yddot_des = A_a * sin(w * t_ref);
ydot_des = (A_a/w) * (1 - cos(w * t_ref));
y_des = -(A_a/w^2) * sin(w*t_ref) + (A_a/w) * t_ref;



x0 = [0; 0; 0; 0];  % Initial state
x0_aug = [x0; 0];

% === REPLACE YOUR FOR LOOP WITH THIS === %
% Standard integral action tracking
closed_loop_dynamics_aug = @(t, x_aug) A_aug*x_aug + B_aug*(-K*x_aug) + [ones(4,1)*0; 1] * interp1(t_ref, y_des, t);
[t_sim, X_aug] = ode45(closed_loop_dynamics_aug, t_ref, x0_aug);  % SIMULATE!

% Extract output and compute control history
X = X_aug(:, 1:4)';  % Original states
Y = X(2, :);         % Output x2
integral_states = X_aug(:, 5)';

U = ones(1, length(t_sim)) * 0;
for k = 1:length(t_sim)
    U(k) = -K * X_aug(k, :)';  % Use augmented state directly
end
t = t_sim';  % Use the time points from ODE45
% === END REPLACEMENT === %

inte_U = trapz(U);
disp(inte_U)


figure;
subplot(3,1,1);
plot(t, y_des, 'r--', t, Y, 'b');
legend('desired', 'actual');
ylabel('Position');

subplot(3,1,2);
plot(t, U);
ylabel('control input');

subplot(3,1,3);
plot(t, y_des - Y);
ylabel('tracking error')