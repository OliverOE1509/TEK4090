% After computing theoretical and sample variances
theoretical_var_WLS = 1 / (N * (3.333e-5 / 1^2) + N * (3.333e-5 / 16^2));
theoretical_var_LS = ((1^2 + 16^2)/2) / (2 * N * 3.333e-5);

sample_var_WLS = var(R_WLS);
sample_var_LS = var(R_LS);

disp(theoretical_var_WLS)
disp(sample_var_WLS)

% Since variances can differ by orders of magnitude
figure;
semilogy(1, theoretical_var_WLS, 'bo', 'MarkerSize', 10, 'LineWidth', 2);
hold on;
semilogy(1, sample_var_WLS, 'bx', 'MarkerSize', 10, 'LineWidth', 2);
semilogy(2, theoretical_var_LS, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
semilogy(2, sample_var_LS, 'rx', 'MarkerSize', 10, 'LineWidth', 2);

xlim([0.5, 2.5]);
set(gca, 'XTick', [1, 2], 'XTickLabel', {'WLS', 'LS'});
ylabel('Variance (log scale)');
title('Variance Comparison (Logarithmic Scale)');
legend('Theoretical WLS', 'Empirical WLS', 'Theoretical LS', 'Empirical LS', ...
       'Location', 'northwest');
grid on;