clear; clc; close all;

h = 0.5; % bar length
M = 0.1; % propeller mass
m = 0.4; % bar mass
b = 1; % viscosity coefficient
Ts = 0.1; % sampling time
g = 9.8; % gravity acceleration
Jeq = (M + m/3)*h^2; % Resultant inertia momentum
Par = [h M m b g Jeq];

% Simulation
Tsim = 10;
num = floor(Tsim/Ts);

x0 = [0*pi/180 0]'; % initial conditions
tspan = Ts:Ts:Tsim; % integration interval
u = 1*ones(1,length(tspan)); % magnitude of step input
options = odeset('RelTol',1e-12,'AbsTol',1e-12*ones(1,2)); % making the solutions very accurate
[t,x] = ode45(@(t,x) propeller(t,x,Par,tspan,u),tspan,x0,options); % solves the differential equation system
theta = x(:,1); % takes the first column (theta)

subplot(2,1,1)
plot(t,theta*180/pi,"LineWidth",1)
ylabel("Agular Position")
grid on

subplot(2,1,2)
plot(t,u,"LineWidth",1)
ylabel("Control Signal")
xlabel("Time")
grid on