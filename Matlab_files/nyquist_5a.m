s = tf('s');
P = (s^2 + 1)/((s^2 + 2)*(s - 2));

% Frequency range, avoiding the pole at sqrt(2)
w1 = logspace(-2, log10(sqrt(2))-0.05, 400);
w2 = logspace(log10(sqrt(2))+0.05, 2, 400);

figure;
nyquist(P, w1); hold on;
nyquist(P, w2);

grid on;
title('Nyquist plot of P(s) = (s^2+1)/((s^2+2)(s-2))');