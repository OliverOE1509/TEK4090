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

disp(zeros)

T = 10;
w = 2 * pi / T;
A_a = w / T;
t = 0:0.01:T;


y_des_func = @(t) -(A_a/w^2) * sin(w*t) + (A_a/w) * t;
yd_des_func = @(t) (A_a/w) * (1 - cos(w*t));
ydd_des_func = @(t) A_a * sin(w*t);

% Calculate the desired trajectories at all time points
y_des = y_des_func(t);
ydot_des = yd_des_func(t);
yddot_des = ydd_des_func(t);

disp(y_des)

return


internal_dynamics = @(t, x1_vec) [x1_vec(2); 0.1 * x1_vec(2) + x1_vec(1) - 12 * ydd_des_func(t) - 0.1 * yd_des_func(t) - y_des_func(t)];
x1_initial = [0;0];
[t_internal, X1_internal] = ode45(internal_dynamics, t, x1_initial);


x1_ref = X1_internal(:, 1)';
x1dot_ref = X1_internal(:, 2)';
 

x1ddot_ref = gradient(x1dot_ref, t);

% Calculate inverse input
u_inv = 12*x1ddot_ref + yddot_des + 0.1*x1dot_ref - 0.1*ydot_des + x1_ref - y_des;

% ===== COMPLETE REFERENCE STATE =====
x_ref = [x1_ref; y_des; x1dot_ref; ydot_des];

% Plot the internal dynamics
figure;
subplot(2,1,1);
plot(t, x1_ref, 'b', t, y_des, 'r--');
legend('x1_{ref} (internal)', 'x2_{ref} (desired output)');
ylabel('Position');
title('Internal Dynamics - Reference States');

subplot(2,1,2);
plot(t, u_inv);
ylabel('Inverse Input u_{inv}');
xlabel('Time (s)');
title('Feedforward Control Input');

% Display some statistics
fprintf('Max inverse input: %.4f\n', max(abs(u_inv)));
fprintf('Max x1_ref: %.4f\n', max(abs(x1_ref)));

