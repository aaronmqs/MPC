% Author: Aaron Marques
% Date: 04/23/2022
% Adaptive Predictive Control - CARIMA - Restrictions
%%
clear, clc, close all;
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot,'defaultTextInterpreter','latex');

s = tf('s');

%%
% Model information
G1 = 1/(2*s+1);
G2 = -0.5/(3*s+1);
Ts = 0.1;
z = tf('z',Ts);
G1d = c2d(G1,Ts);
G2d = c2d(G2,Ts);
[b1b,a1]  = tfdata(G1d,'v');
[b2b,a2]  = tfdata(G2d,'v');
delta = [1 -1];
alfa = 0; % parâmetro de ajuste
c = conv([1 -alfa],[1 -alfa]);
%c = 1;
a_til = conv(delta,conv(a1,a2));
d = a_til;
b1 = conv(b1b,a2);
b2 = conv(b2b,a1);
[a_til,b1] = eqtflength(a_til,b1);
[a_til,b2] = eqtflength(a_til,b2);
[c,d] = eqtflength(c,d); % necessário para usar o tf2ss

%%
% Space-State Model
[~,~,B1,~] = tf2ss(b1,a_til);
B1= B1';

% Space-State Model
[A,H,B2,~] = tf2ss(b2,a_til);
A = A';
B2 = B2';
H = H';

B = [B1 B2];

[~,~,D,~] = tf2ss(c,d);
D = D';

% Parâmetros de controle
N = 10; % Horizonte de predição
Nu = 1; % Horizonte de controle

%%
% Matrix definitions
F = cell(N,1);
E = cell(N,1);
G = cell(N,Nu);

for i = 1:N
    F{i,1} = H*A^i;
    E{i,1} = H*(A^(i-1))*D;

    for j = 1:Nu
       if i<j
           G{i,j} = zeros(size(H*B));
       else
           G{i,j} = H*(A^(i-j))*B;
       end
    end

end

F = cell2mat(F);
E = cell2mat(E);
G = cell2mat(G);
Q_lamb = eye(size(G,2)); % peso do erro de seguimento à referência no custo
Aux = (G'*G + Q_lamb)\G';
K = Aux(1:2,:);

%%
%Quadprog parameters
Hb = 2*(G'*G + Q_lamb);
Bo = cell(Nu,1);
M = cell(Nu);

for i = 1:Nu

    if i == 1

        Bo{i,1} = eye(size(B,2));

    else

        Bo{i,1} = zeros(size(B,2));

    end

    for j = 1:Nu

        if i == j

            M{i,j} = eye(size(B,2));

        elseif i == j+1

            M{i,j} = -eye(size(B,2));

        else

            M{i,j} = zeros(size(B,2));

        end
    end

end

M = cell2mat(M);
T = inv(M);
Bo = cell2mat(Bo);
TBo = M\Bo;
Ab = [T;-T];
um = 0;
uM = 3;
Umin = um*ones(size(M,1),1);
Umax = uM*ones(size(M,1),1);

%%
% Simulation
Tsim = 20;
Num = floor(Tsim/Ts); % Number of samples
t = Ts:Ts:Num*Ts;

x = cell(1,Num);
y1 = cell(1,Num);
y2 = cell(1,Num);
y = cell(1,Num);
u = cell(1,Num-1);
e = cell(1,Num-1);

% Initial conditions
x{1} = zeros(size(A,2),1);
[y1{1},z1] = filter(b1b,a1,0);
[y2{1},z2] = filter(b2b,a2,0);
% y1{1} = 0;
% y2{1} = 0;
y{1} = y1{1} + y2{1};

for i = 1:Num-1

    e{i} = y{i} - H*x{i};
    f = F*x{i} + E*e{i};
    b = 2*G'*(f - 1);
    
    if i == 1
        TBo_u = zeros(size(Umax));
    else
        TBo_u = TBo*u{i-1};
    end

    Bb = [Umax-TBo_u;TBo_u-Umin];
    Delta_U = quadprog(Hb,b,Ab,Bb); % control effort: Delta_U = [Delta_u(t) Delta_u(t+1) ... Delta_u(t + N - 1)]'
    U = M\Delta_U + TBo_u; % control signal: U = [u(t) u(t+1) ... u(t + N - 1)]'
    u{i} = [U(1) U(2)]'; % u(t) = [u1(t) u2(t)]'

    % process Input Disturbance
    if i > 200
        Up(1) = U(1) + 1;
        Up(2) = U(2) + 1;
    else 
        Up(1) = U(1);
        Up(2) = U(2);
    end

%     It's not possible using LSIM in this case since it'd be needed using past states for each system.    
%     output = lsim(G1d,U(1:2:end)) + lsim(G2d,U(2:2:end));
%     y{i+1} = output(2);

    %y1{i+1} = 0.9512*y1{i} + 0.04877*Up(1);
    [y1{i+1},z1] = filter(b1b(2:end),a1,Up(1),z1);
    %y2{i+1} = 0.9672*y2{i} - 0.01639*Up(2);  
    [y2{i+1},z2] = filter(b2b(2:end),a2,Up(2),z2);

    y{i+1} = y1{i+1} + y2{i+1};    

    x{i+1} = A*x{i} + B*Delta_U(1:size(B,2)) + D*e{i};

end

y = cell2mat(y);
u = cell2mat(u);

subplot(2,1,1)
plot(t,y,'k','LineWidth',3)
ylabel('Output Signal')
legend('y(t)')
grid on
subplot(2,1,2)
stairs(t(1:end-1),u(1,:),'LineWidth',3)
hold on
stairs(t(1:end-1),u(2,:),'LineWidth',3)
ylabel('Control Signals')
xlabel('Time')
legend('u_1(t)','u_2(t)')
grid on
sgtitle("Outpupt and Control Signals for N = " + N + ", Umin = " + um + ", Umax = " + uM)





















