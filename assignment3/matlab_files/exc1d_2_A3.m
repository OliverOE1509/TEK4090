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
R = 0.1;
K = lqr(A, B, Q, R);

yddot_des = A_a * sin(w * t_ref);
ydot_des = (A_a/w) * (1 - cos(w * t_ref));
y_des = -(A_a/w^2) * sin(w*t_ref) + (A_a/w) * t_ref;



x0 = [0; 0; 0; 0];  % Initial state

% === REPLACE YOUR FOR LOOP WITH THIS === %
closed_loop_dynamics = @(t, x) A*x + B*(-K*(x - [0; interp1(t_ref, y_des, t); 0; interp1(t_ref, ydot_des, t)]));
[t_sim, X] = ode45(closed_loop_dynamics, t_ref, x0);


% Extract output and compute control history
Y = X(:, 2)';


fprintf('Length of t_sim')
disp(length(t_sim))



%U = zeros(1, length(t_sim));  % Use number of rows instead of length

%U = zeros(1, 1001);
U = ones(1, length(t_sim)) * 0;
for k = 1:length(t_sim)
    x_ref = [0; interp1(t_ref, y_des, t_sim(k)); 0; interp1(t_ref, ydot_des, t_sim(k))];
    U(k) = -K * (X(k, :)' - x_ref);
end
t = t_sim';  % Use the time points from ODE45
% === END REPLACEMENT === %

y_error = y_des - Y;
inte_y_error = trapz(y_error);
inte_U = trapz(U);
fprintf('U control effort')
disp(inte_U)
fprintf('Y error sum ')
disp(inte_y_error)

figure;
subplot(3,1,1);
plot(t, interp1(t_ref, y_des, t), 'r--', t, Y, 'b');  % Interpolate desired for plotting
legend('desired', 'actual');
ylabel('Position');
title('LQR Tracking with ODE45 Integration');

subplot(3,1,2);
plot(t, U);
ylabel('control input');

subplot(3,1,3);
plot(t, y_des - Y);
ylabel('tracking error')