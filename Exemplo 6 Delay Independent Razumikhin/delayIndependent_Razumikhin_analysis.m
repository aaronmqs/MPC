clear
clc;
close all;

% dx/dt(t) = A0x(t) + Bu(t-h)
% u(t) = Kx(t)
% G(s) = 2/(s^2 + 3s + 2)
% The static gain is 1, thus K = [1 0] makes the system marginally stable
% That's because the maximum magnitude of the Bode diagram is 0dB (1)
A0 = [-3 1;-2 0];
B = [0 2]';
K = -[0.9 0]; % The system is stable
% K = -[1.1 0]; % The system is unstable
A1 = B*K;

P = sdpvar(size(A0,1));
S = sdpvar(size(A0,1));

LMI = [ P >= 0;
        S >= 0;
        [P*A0+A0'*P+S,P*A1;
        A1'*P,-S] <= 0];

optimize(LMI);
checkset(LMI)
Po = value(P);
So = value(S);

sys = ss(A0,B,K,0);

subplot(1,2,1)
bode(sys)
grid on
subplot(1,2,2)
step(sys)
grid on