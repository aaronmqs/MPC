clear
clc;
close all;

A = [-3 1;-2 0];
B = [0 2]';

Q = sdpvar(2, 2, 'symmetric', 'real');
Y = sdpvar(1, 2,'full','real');
W = sdpvar(2, 2, 'full', 'real');

LMI1 = (Q >= 0);

LMI3 = (W >= 0);

LMI2 = [(A*Q + Q*A' + W), B*Y;
    Y'*B', -W] <= 0;

LMIs = [LMI1, LMI2; LMI3];

options = sdpsettings('solver','sedumi','verbose',0);
optimize(LMIs);

Q = value(Q);
Y = value(Y);
W = value(W);

P = inv(Q);
alpha = W/P
K = Y*P

sys = ss(A,B,K,0);
bode(sys)






























% Gabriel's code
% 
% Q = sdpvar(2, 2, 'symmetric', 'real');
% Y = sdpvar(1, 2,'full','real');
% W = sdpvar(2, 2, 'full', 'real');
% 
% LMI1 = (Q >= 0);
% 
% LMI3 = (W >= 0);
% 
% LMI2 = [(A*Q + Q*A' + W), B*Y;
%     Y'*B', -W] <= 0;
% 
% LMIs = [LMI1, LMI2; LMI3];
% 
% options = sdpsettings('solver','sedumi','verbose',0);
% optimize(LMIs, W);
% 
% Q = value(Q);
% Y = value(Y);
% W = value(W);
% 
% P = inv(Q);
% alpha = W/P
% K = Y*P