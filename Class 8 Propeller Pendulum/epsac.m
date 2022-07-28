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
N = 10; % prediction horizon
alpha = 0.4;
delta = [1 -1];
c = [1 -alpha];

%% Simulation
Tsim = 3;
num = ceil(Tsim/Ts);
%N = num;
t = Ts:Ts:N*Ts;

% Reference
ref = (50*pi/180)*ones(1,num);
% ref = [];
% for i = 1:18
%     ref = [ref, 10*i*(pi/180)*ones(1,50)]; 
% end

% Equilibrium point
theta_eq = 10*pi/180; % desired equilibrium angle (intial position)
ueq = M*g*sin(theta_eq)*h + m*g*sin(theta_eq)*h/2; % input needed to keep the initial position

% initial conditions
y1 = 40*pi/180; 
y2 = 40*pi/180; 
y(1) = 40*pi/180; % Measured output
u1 = 0;
u2 = 0;
[e1,ze1] = filter(delta,c,0);
[np,znp] = filter((1-alpha),[1 -1],0);

for k = 1:num

    x(k) = 2*y1 - y2 + (Ts^2/Jeq)*(-M*g*h*sin(y2) - m*g*h*sin(y2)/2 - (b/Ts)*(y1 - y2) + u2);
    n(k) = y(k) - x(k);
    [e(k),ze1] = filter(delta,c,n(k),ze1);

    [np,znp] = filter((1-alpha),[1 -1],e(k),znp);
    n = np*ones(1,N);

    x11 = y(k);
    x12 = y1;
    x21 = y(k);
    x22 = y1;
    a = 0.0001;
    
    for i =1:N
        
        % free response
        x1(i) = 2*x11 - x12 + (Ts^2/Jeq)*(-M*g*h*sin(x12) - m*g*h*sin(x12)/2 - (b/Ts)*(x11 - x12) + u1);
        y11(i) = x1(i) + n(i);
        x12 = x11;
        x11 = y11(i);
        
        % step response
        x2(i) = 2*x21 - x22 + (Ts^2/Jeq)*(-M*g*h*sin(x22) - m*g*h*sin(x22)/2 - (b/Ts)*(x21 - x22) + u1 + a);
        y22(i) = x2(i) + n(i);
        x22 = x21;
        x21 = y22(i);

    end

    G = (y22 - y11)'/a; % unit step response
    uo = (G'*G)\G'*(ref(k)*ones(N,1) - y11'); % optimal input control signal
    u(k) = uo + u1; % control signal
    
    % Saturation
    if u(k) > 2.5
        u(k) = 2.5;
    end

    y2 = y1;
    y1 = y(k);
    
    %% Reading the next measure of the real process
    % Uncertainty
    if k >= 100

        y(k+1) = 2*y1 - y2 + (Ts^2/Jeq)*(-M*1.5*g*h*sin(y2) - m*g*h*sin(y2)/2 - (b/Ts)*(y1 - y2) + u1 - 0.5);
    else
        y(k+1) = 2*y1 - y2 + (Ts^2/Jeq)*(-M*1.5*g*h*sin(y2) - m*g*h*sin(y2)/2 - (b/Ts)*(y1 - y2) + u1);
    end
    
    % noise
    if k>= 200
        y(k+1) = y(k+1) + 10^-3*rand;
    end
    %%
    
    u2 = u1;
    u1 = u(k);
    
end

subplot(2,1,1)
plot(y*180/pi)
hold on
plot(ref*180/pi)
grid on

subplot(2,1,2)
plot(u)
grid on

