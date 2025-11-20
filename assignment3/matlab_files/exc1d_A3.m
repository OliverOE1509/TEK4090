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

T = 100;
w = 2 * pi / T;
A_a = w / T;
t = 0:0.01:T;



Q = diag([1, 5000, 1, 50]);
R = 0.01;
K = lqr(A, B, Q, R);

yddot_des = A_a * sin(w * t);
ydot_des = (A_a/w) * (1 - cos(w * t));
y_des = -(A_a/w^2) * sin(w*t) + (A_a/w) * t;



x0 = [0; 0; 0; 0];  % Initial state
x = x0;
X = []; U = []; Y = [];

for k=1:length(t)
    x_ref = [0; y_des(k); 0; ydot_des(k)];
    u_fb = -K * (x - x_ref);
    
    X = [X, x];
    U = [U, u_fb];
    Y = [Y, x(2)];

    if k < length(t)
        dt = t(k+1) - t(k);
        xdot = A*x + B*u_fb;
        x = x + xdot * dt;
    end
end

y_error = y_des - Y;
inte_y_error = trapz(y_error);
inte_U = trapz(U);
fprintf('U control effort')
disp(inte_U)
fprintf('Y error sum ')
disp(inte_y_error)




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