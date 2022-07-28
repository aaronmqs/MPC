clear; clc; close all;

Tsim = 29.98;

%% Plant model
h = 0.5; % bar length
M = 0.1; % propeller mass
m = 0.4; % bar mass
b = 0.1; % viscosity coefficient
k = 1; % spring constant
Ts = 0.01; % sampling time
g = 9.81; % gravity acceleration
Jeq = (M + m/3)*h^2; % Resultant inertia momentum

%% OpenSim Data

u = load("C:\Users\aaron\OneDrive\Área de Trabalho\Backup do Drive\2022.1\Courses\MPC\Codes\Final Project\opensim data\muscleexcitation.mat");
x = load("C:\Users\aaron\OneDrive\Área de Trabalho\Backup do Drive\2022.1\Courses\MPC\Codes\Final Project\opensim data\elbowangle.mat");
t = load("C:\Users\aaron\OneDrive\Área de Trabalho\Backup do Drive\2022.1\Courses\MPC\Codes\Final Project\opensim data\time.mat");

u = u.excitationplot-0.01;
x = x.elbowangleplot;
x = x - x(1);
% t = t.controlTimeplot;
t = 0.02*(0:1499);

%% Simulink Data

time = t';
input_data = u';
input = [time,input_data];

%% Plant model plot

simulation = sim("SystemSimulation.slx",Tsim);
tout = simulation.tout;
theta = simulation.theta;

subplot(311)
plot(t,u)
title('Muscle excitation')
grid on
subplot(312)
plot(tout,theta)
title('Elbow angle - Phenomenologiacal model')
grid on

%% Optimized plant model
h = 0.3632;
M = 0.0857;
m = 0.2853;
b = 0.0686;
k = -0.0056;
g = 0.8234;
Jeq = 5.4579e-4;

simulation = sim("SystemSimulation.slx",Tsim);
tout = simulation.tout;
theta = simulation.theta;

%% Expeerimental data

subplot(313)
plot(t,x)
hold on
plot(tout,theta)
legend("Experimental data","Optimized model")
title('Elbow angle - Compared results')
grid on














