s = tf('s');
P = (s^2+1)/(s^3+2*s^2+2*s+1);
C = 1;
pc = P * C;

% Use nyquistplot with specified frequency range
figure;
nyquistplot(pc, {0.01, 1000}); % Adjust frequency range as needed
grid on;
title('Nyquist Plot of PC(s) = (s^2+1)/(s^3+2*s^2+2*s+1)');