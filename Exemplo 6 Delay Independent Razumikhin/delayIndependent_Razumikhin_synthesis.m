clear
clc;
close all;

% dx/dt(t) = A0x(t) + Bu(t-h)
% u(t) = Kx(t)
% G(s) = 2/(s^2 + 3s + 2)
% The static gain is 1, thus K = [1 0] makes the system marginally stable
% That's because the maximum magnitude of the Bode diagram is 0dB (1)
A0 = [-3 1;-2 0];
A1 = [0 2]';
C = eye(2);
D = 0;

% Q = alpha/P
% Y = K/P
P = sdpvar(size(A0,1),size(A0,1),'symmetric');
Q = sdpvar(size(A0,1),size(A0,1),'symmetric');
Y = sdpvar(size(A1,2),size(A0,2));

LMI = [ P >= 0;
        Q >= 0;
        [P*A0+A0'*P+Q, A1*Y;
        Y'*A1', -Q] <= 0;
    ];

optimize(LMI);
checkset(LMI);

P = value(P);
Q = value(Q);
Y = value(Y);

K = Y*P;

sys_OL = ss(A0,A1,C,0);
sys_CL = ss(A0+A1*K,A1,K,0);

subplot(2,1,1)
bode(sys_OL)
grid on

subplot(2,1,2)
bode(sys_CL)
grid on
