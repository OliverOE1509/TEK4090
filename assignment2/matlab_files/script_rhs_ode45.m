clear all; close all; clc;

% simulation timespan
dt=0.1;
t_space = [0:dt:10];

% initial condition
x_0 = [10;10];

% sinusoidal input with zero-order-hold
u = cos(t_space);
x_array = zeros(2,length(t_space));

x_array(:,1) = x_0;

for ii = 1:length(t_space)-1

    t_int = [t_space(ii), t_space(ii+1)];

    %simulate ODE with input
    [t_out,y_out] = ode45(@fn_rhs,t_int,x_array(:,ii),[],u(ii));
    x_array(:,ii+1) = y_out(end,:)';

end

 plot(t_space, x_array(1,:)), hold on
 plot(t_space, x_array(2,:))