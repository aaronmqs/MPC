clear; clc; close all;


%% Plant model
h = 0.5; % bar length
M = 0.1; % propeller mass
m = 0.4; % bar mass
b = 0.1; % viscosity coefficient
Ts = 0.01; % sampling time
g = 9.81; % gravity acceleration
Jeq = (M + m/3)*h^2; % Resultant inertia momentum

%% Control parameters
N = 1000; % prediction horizon

%% Simulation
Tsim = 10;
num = ceil(Tsim/Ts);
%N = num;
t = Ts:Ts:N*Ts;

a = 0.01; % input magnitude
ub = a*ones(1,N); % base input

% Equilibrium point
theta_eq = 15*pi/180; % desired equilibrium angle (intial position)
ueq = M*g*sin(theta_eq)*h + m*g*sin(theta_eq)*h/2; % input needed to keep the initial position

u = ub + ueq*ones(1,N);

% initial conditions
theta1 = theta_eq;
theta2 = theta_eq;

for k = 1:N
    theta(k) = 2*theta2 - theta1 + (Ts^2/Jeq)*(-M*g*h*sin(theta1) - m*g*h*sin(theta1)/2 - (b/Ts)*(theta2 - theta1) + u(k));
    theta1 = theta2;
    theta2 = theta(k);
    % theta(k+2) = 2*theta(k+1) - theta(k) + (Ts^2/Jeq)*(-M*g*h*sin(theta(k)) - m*g*h*sin(theta(k))/2 - (b/Ts)*(theta(k+1) - theta(k)) + u(k));
end

subplot(3,1,1)
plot(Ts*(1:N),theta*180/pi,"LineWidth",1)
title("Ts = " + Ts + " and N = " + N + " ")
ylabel("Agular Position (degrees)")
grid on

subplot(3,1,2)
stairs(Ts*(1:N),u,"LineWidth",1)
title("a = " + a + "")
ylabel("Control Signal (N.m)")
grid on

subplot(3,1,3)
plot(Ts*(1:N),((theta - theta_eq)*180/pi)/a,"LineWidth",1)
title("Step Response")
ylabel("Agular Position (degrees)")
xlabel("Time (s)")
grid on

