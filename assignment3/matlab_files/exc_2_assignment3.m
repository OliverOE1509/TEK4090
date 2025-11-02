sigma_good_voltmeter = 1;
sigma_bad_voltmeter = 16;
N = 100;
R0 =  1250; % Ohms
i_uniform = unifrnd(-0.01, 0.01, N, 1);
nu1 = normrnd(0, sigma_good_voltmeter, N, 1);
nu2 = normrnd(0, sigma_bad_voltmeter, N, 1);

voltage_good = R0 * (i_uniform + nu1);
voltage_bad = R0 * (i_uniform + nu2);

w1 = var(voltage_good);
w2 = var(voltage_bad);
voltage_cat = horzcat(voltage_good', voltage_bad');

R_hat_WLS = 0;
for t = 1;2*N
    R_hat_WLS += ( (i_uniform * voltage_cat(t) ) /   

