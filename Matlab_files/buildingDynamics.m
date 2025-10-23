function dzdt = buildingDynamics(t, z, A, B, K, x_bar, b_ambient, D_disturbance, x_a, d_func)
    % Control law: u = -K*z
    u = -K * z;
    
    % CORRECT: System dynamics in z-coordinates
    % dz/dt = A*(z + x_bar) + B*u + b_ambient*x_a(t) + D_disturbance*d_func(t) - A*x_bar
    % But A*x_bar is constant offset, so we keep A*z
    dzdt = A*z + B*u + b_ambient*(x_a(t) - 20) + D_disturbance*d_func(t);
end