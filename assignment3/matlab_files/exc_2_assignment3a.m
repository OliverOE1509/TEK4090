sigma_good_voltmeter = 1;
sigma_bad_voltmeter = 16;
N = 100;
R0 =  1250; % Ohms

function R_hat_WLS = generate_Rhat_WLS(N, R0, sigma_good, sigma_bad)
    i0 = -0.01 + 0.02 * rand(N, 1);    
    nu1 = sigma_good * randn(N, 1);
    nu2 = sigma_bad* randn(N, 1);
 
    voltage_good = R0 * (i_uniform + nu1); % This is the u
    voltage_bad = R0 * (i_uniform + nu2);    

    u0 = R0 * i0;
    u_good = u0(:N) + randn(N, 1); % Sigma 2 = 1
    u_bad = u0(N+1: 2 * N) + 4 * randn(N, 1); %sigma"2 = 16
    
    w1 = var(voltage_good);
    w2 = var(voltage_bad);
    voltage_cat = horzcat(voltage_good', voltage_bad');
    numerator = 0;
    denominator = 0;
    
    for t = 1:(2*N)
        if t <= N
            % First N points: good voltmeter measurements
            current_t = i_uniform(t);
            voltage_t = voltage_cat(t);
            weight = w1;
        else
            % Last N points: bad voltmeter measurements  
            current_t = i_uniform(t - N);  % Adjust index for i_uniform
            voltage_t = voltage_cat(t);
            weight = w2;
        end
        
        numerator = numerator + (current_t * voltage_t) / weight;
        denominator = denominator + (current_t^2) / weight;
    end
    
    R_hat_WLS = numerator / denominator;
end

function R_hat_LS = generate_Rhat_LS(N, R0, sigma_good, sigma_bad)
    i_uniform = -0.01 + 0.02 * rand(N, 1);    
    nu1 = sigma_good * randn(N, 1);
    nu2 = sigma_bad* randn(N, 1);
 
    voltage_good = R0 * (i_uniform + nu1); % This is the u
    voltage_bad = R0 * (i_uniform + nu2);    
    
    voltage_cat = horzcat(voltage_good', voltage_bad');
    numerator = 0;
    denominator = 0;
    
    for t = 1:(2*N)
        if t <= N
            % First N points: good voltmeter measurements
            current_t = i_uniform(t);
            voltage_t = voltage_cat(t);
        else
            % Last N points: bad voltmeter measurements  
            current_t = i_uniform(t - N);  % Adjust index for i_uniform
            voltage_t = voltage_cat(t);
        end
        
        numerator = numerator + (current_t * voltage_t);
        denominator = denominator + (current_t^2);
    end
    
    R_hat_LS = numerator / denominator;
end



R_WLS = zeros(1e5, 1);
R_LS = zeros(1e5, 1);
for i = 1:1e5
    R_WLS(i) = generate_Rhat_WLS(N, R0, sigma_good_voltmeter, sigma_bad_voltmeter);
    R_LS(i) = generate_Rhat_LS(N, R0, sigma_good_voltmeter, sigma_bad_voltmeter);
end


disp(var(R_WLS))
disp(var(R_LS))


return
figure;
h1 = histogram(R_WLS, 'FaceAlpha', 0.5, 'DisplayName', 'WLS');
hold on;
h2 = histogram(R_LS, 'FaceAlpha', 0.5, 'DisplayName', 'LS');
xlabel('Resistance (Ohms)');
ylabel('Frequency');
title('Comparison of R_{LS} and R_{WLS} Estimates');
legend;
grid on;