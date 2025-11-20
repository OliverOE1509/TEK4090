N = 100;
a = 0.2;
b = -0.25;
%eps = normrnd(0, 0.5, N, 1);
eps = randn(N, 1);
eps = 0.5 * eps;

function u = generate_random_u_of_1(N)
    u = 2 * rand(N, 1) - 1; 
end


function [y_, Phi] = generate_y(N, a, b)
    u = generate_random_u_of_1(N);
    eps = normrnd(0, 0.5, N, 1);
    y_ = zeros(N, 1);
    y_(1) = 0;
    Phi = zeros(N-1, 2);
    for i = 2:N
        y_(i) = a * y_(i-1) + b * u(i-1) + eps(i);
        Phi(i-1, 1) = y_(i-1);
        Phi(i-1, 2) = u(i-1);
    end
    y_ = y_(2:end);
end

[y_out, Phi] = generate_y(N, a, b);

function theta = generate_theta(Phi, y)
    theta = (Phi' * Phi) \ (Phi' * y);
end
theta = generate_theta(Phi, y_out);
disp(theta)


% Generate a and b 1e5 times to get a histogram distribution
num_experiments = 1e5;
a_hat = zeros(num_experiments, 1);
b_hat = zeros(num_experiments, 1);

for i = 1:1e5
    [y_out, Phi] = generate_y(N, a, b);
    theta = generate_theta(Phi, y_out);
    a_hat(i) = theta(1);
    b_hat(i) = theta(2);
end



% Plot histograms
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
