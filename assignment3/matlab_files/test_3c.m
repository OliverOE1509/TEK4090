N = 100;
a = 0.2;
b = -0.25;
eps = normrnd(0, 0.5, N, 1);

function [y_, Phi] = generate_y(N, a, b, which_u)
    switch which_u
        case 1
            %u = 2*randi([0, 1], N, 1) - 1;
            u = ones(N, 1);
            %‘u = zeros(N, 1);
        case 2
            u = 2 * rand(N, 1) - 1; 
        case 3
            u = randn(N, 1);
        otherwise
            error("Invalid input method. Use 1, 2 or 3")
    end
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

which_u = 1; % Example input method
[y_out, Phi] = generate_y(N, a, b, which_u);

function theta = generate_theta(Phi, y)
    theta = (Phi' * Phi) \ (Phi' * y);
end
theta = generate_theta(Phi, y_out);

Phiss = transpose(Phi) * Phi;
d = det(Phiss);
disp(Phiss)
disp(d)


