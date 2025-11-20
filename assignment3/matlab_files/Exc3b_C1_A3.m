N = 100;
a = 0.2;
b = -0.25;
%eps = normrnd(0, 0.5, N, 1);
eps = randn(N, 1);
eps = 0.5 * eps;

function [y_, Phi] = generate_y(N, a, b, which_u)
    switch which_u
        case 1
            u = 2*randi([0, 1], N, 1) - 1;
            %u = ones(N, 1);
        case 2
            u = 2 * rand(N, 1) - 1; 
        case 3
            u = randn(N, 1);
        otherwise
            error("Invalid input method. Use 1, 2 or 3")
    end
    eps = normrnd(0, 0.5, N, 1);
    y_ = zeros(N, 1);
    y_(1) = 100;
    Phi = zeros(N-1, 2);
    for i = 2:N
        y_(i) = a * y_(i-1) + b * u(i-1) + eps(i);
        Phi(i-1, 1) = y_(i-1);
        Phi(i-1, 2) = u(i-1);
    end
    y_ = y_(2:end);
end

function theta = generate_theta(Phi, y)
    theta = (Phi' * Phi) \ (Phi' * y);
end
theta = generate_theta(Phi, y_out);


% Generate a and b 1e5 times to get a histogram distribution
num_experiments = 1e5;
a_hat1 = zeros(num_experiments, 1);
b_hat1 = zeros(num_experiments, 1);

a_hat2 = zeros(num_experiments, 1);
b_hat2 = zeros(num_experiments, 1);

a_hat3 = zeros(num_experiments, 1);
b_hat3 = zeros(num_experiments, 1);


d = struct();
d.method1 = [1, 2];
d.method2 = [1, 2];
d.method3 = [1, 2];
methods = fieldnames(d);
disp(methods{1})

% Create one big figure with 3x2 subplots
figure('Position', [100, 100, 1200, 800]);

% Method names for titles
method_names = {'Random +1 -1', 'Uniform [-1,1]', 'Gaussian N(0,1)'};

% Loop through all methods
for method_idx = 1:length(methods)
    % Get current method arrays
    a_hat_array = eval(['a_hat' num2str(method_idx)]);
    b_hat_array = eval(['b_hat' num2str(method_idx)]);
    
    % Generate estimates for this method
    for i = 1:num_experiments
        [y_out, Phi] = generate_y(N, a, b, method_idx);
        theta = generate_theta(Phi, y_out);
        a_hat_array(i) = theta(1);
        b_hat_array(i) = theta(2);
    end
    
    % Store back in original arrays
    eval(['a_hat' num2str(method_idx) ' = a_hat_array;']);
    eval(['b_hat' num2str(method_idx) ' = b_hat_array;']);
    
    % Plot a histogram (left column)
    subplot(3, 2, 2*method_idx - 1);
    histogram(a_hat_array, 50, 'FaceColor', 'blue', 'FaceAlpha', 0.7);
    hold on;
    xline(a, 'r--', 'LineWidth', 2, 'Label', 'True a');
    title(['hat(a) - ' method_names{method_idx}]);
    xlabel('Estimated a');
    ylabel('Frequency');
    xlim([-0.1, 0.5]);
    grid on;
    legend('Estimates', 'True value');
    
    % Plot b histogram (right column)
    subplot(3, 2, 2*method_idx);
    histogram(b_hat_array, 50, 'FaceColor', 'green', 'FaceAlpha', 0.7);
    hold on;
    xline(b, 'r--', 'LineWidth', 2, 'Label', 'True b');
    title(['hat(b) - ' method_names{method_idx}]);
    xlabel('Estimated b');
    ylabel('Frequency');
    xlim([-0.5, 0]);
    grid on;
    legend('Estimates', 'True value');
end

variances_a = zeros(3,1);
variances_b = zeros(3,1);

for method_idx = 1:length(methods)
    a_hat_array = eval(['a_hat' num2str(method_idx)]);
    b_hat_array = eval(['b_hat' num2str(method_idx)]);
    
    var_a = var(a_hat_array);
    var_b = var(b_hat_array);
    variances_a(method_idx) = var_a;
    variances_b(method_idx) = var_b;
    
    fprintf('%-15s %-12.6f %-12.6f\n', method_names{method_idx}, var_a, var_b);
end

% Find best methods
[~, best_a] = min(variances_a);
[~, best_b] = min(variances_b);

fprintf('\n=== RESULTS ===\n');
fprintf('Lowest variance for a: %s\n', method_names{best_a});
fprintf('Lowest variance for b: %s\n', method_names{best_b});



%{

for i = 1:1e5
    [y_out, Phi] = generate_y(N, a, b);
    theta = generate_theta(Phi, y_out);
    a_hat(i) = theta(1);
    b_hat(i) = theta(2);
end

figure;
for i = 1:3
    subplot(1,2,);
    histogram(a_hat, 50);
    title('Histogram of \hat{a} estimates');
    xlabel('Estimated a');
    ylabel('Frequency');
    hold on;
    xline(a, 'r--', 'LineWidth', 2, 'Label', 'True a');
    legend('Estimates', 'True value');
    
    subplot(1,2,2);
    histogram(b_hat, 50);
    title('Histogram of \hat{b} estimates');
    xlabel('Estimated b');
    ylabel('Frequency');
    hold on;
    xline(b, 'r--', 'LineWidth', 2, 'Label', 'True b');
    legend('Estimates', 'True value');

figure;

subplot(1,2,1);
histogram(a_hat, 50);
title('Histogram of \hat{a} estimates');
xlabel('Estimated a');
ylabel('Frequency');
hold on;
xline(a, 'r--', 'LineWidth', 2, 'Label', 'True a');
legend('Estimates', 'True value');

subplot(1,2,2);
histogram(b_hat, 50);
title('Histogram of \hat{b} estimates');
xlabel('Estimated b');
ylabel('Frequency');
hold on;
xline(b, 'r--', 'LineWidth', 2, 'Label', 'True b');
legend('Estimates', 'True value');
%}