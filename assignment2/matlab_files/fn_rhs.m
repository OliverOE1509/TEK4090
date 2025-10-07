function [x_dot] = fn_rhs(t, x, u)
%fn_rhs: a `right-hand-side function' of the ode x_dot = f(t,x,u)

x_dot = exp(-t)*[cos(x(1)); sin(x(2))] + exp(-0.1*t)*[0;1]*u;

end