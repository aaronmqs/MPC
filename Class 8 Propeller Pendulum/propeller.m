function dx = propeller(t,x,Par,tspan,u)

h = Par(1);
M = Par(2);
m = Par(3);
b = Par(4);
g = Par(5);
Jeq = Par(6);
u = interp1(tspan,u,t);

% x(1) = theta
% x(2) = theta dot

dx = [
    x(2);
    (-M*g*h*sin(x(1)) -m*g*h*sin(x(1))/2 - b*x(2) + u)/Jeq;
    ];