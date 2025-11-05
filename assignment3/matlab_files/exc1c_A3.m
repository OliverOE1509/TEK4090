T = 10;
w = 2 * pi / T;
A = w / T;
t = 0:0.01:T;

yddot_des = A * sin(w * t);
ydot_des = (A/w) * (1 - cos(w * t));
y_des = -(A/(w^2)) * sin(w*t) + (A/w) * t;

figure;
subplot(3,1,1);
plot(t, y_des);
ylabel('position y(t)');
grid on;

subplot(3,1,2);
plot(t, ydot_des);
ylabel("velocity y'(t)");
grid on;

subplot(3,1,3);
plot(t, yddot_des);
ylabel("acceleration y''(t)")
grid on;
