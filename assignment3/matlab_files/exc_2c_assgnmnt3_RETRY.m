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
    denominator = 0;
    
    for t = 1:(2*N)
        if t <= N
            % First N points: good voltmeter measurements
            current_t = i_uniform(t);
            weight = w1;
        else
            % Last N points: bad voltmeter measurements  
            current_t = i_uniform(t - N);  % Adjust index for i_uniform
            weight = w2;
        end
        
        numerator = 1;
        denominator = denominator + (current_t^2) / weight;
    end
    
    R_hat_WLS = numerator / denominator;
end

