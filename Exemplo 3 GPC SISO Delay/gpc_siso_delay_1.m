clc, clear, close all;

% Parâmetros do processo
tau1 = 20;
tau2 = 50;
Ts = 1; % tempo de amostragem
z = tf('z',Ts);
s = tf('s');
L = 10;
Ld = ceil(L/Ts);
P = exp(-L*s)/((tau1*s + 1)*(tau2*s + 1));
Pd  = c2d(P,Ts);
[num,den]  = tfdata(Pd,'v');

% Parâmetros do modelo CARIMA
a = den;
b = num;
alfa = 0.9; % parâmetro de ajuste
c = conv([1 -alfa],[1 -alfa]);
delta = [1 -1];
a_til = conv(delta,a);
d = a_til;
[a_til,b] = eqtflength(a_til,b);
[c,d] = eqtflength(c,d); % necessário para usar o tf2ss

% Space-State Model
[A,H,B,~] = tf2ss(b,a_til);
A = A';
B = B';
H = H';

[~,~,D,~] = tf2ss(c,d);
D = D';

% Parâmetros de controle
N = 1; % Horizonte de predição
Nu = 1; % Horizonte de controle

% FIR Filter
FFilter = cell(1,Ld+1);
for i = 1:Ld+1
    if i == 1
        FFilter{i} = zeros(size(A*B));
    else
        FFilter{i} = A^(i-2)*B; % *z^-(i-1)
    end
end

FFilter = cell2mat(FFilter);

% Matrix definitions
F = cell(N,1);
G = cell(N,Nu);

for i = 1:N
    F{i,1} = H*A^i;

    for j = 1:Nu
       if i<j
           G{i,j} = zeros(size(H*B));
       else
           G{i,j} = H*(A^(i-j))*B;
       end
    end

end

Ftil = (eye(size(A)) - (z^-(Ld))*A^(Ld))*(z*eye(size(A)) - A)\B;

F = cell2mat(F);
G = cell2mat(G);
lambda = 1;
delta = 1;
Qlambda = lambda*eye(size(G,2)); % peso do esforço de controle no custo
Qdelta = delta*eye(size(G,1)); % peso do erro de seguimento à referência no custo 

Aux = (G'*Qdelta*G + Qlambda)\G'*Qdelta;
K = Aux(size(B,2),:);
Kr = sum(K);

% simulação no simulink

Tsim = 800;
sim1 = sim('gpc_siso_delay_1_sim');
y = sim1.y;
u = sim1.u;
t = 0:Ts:Tsim;

figure
step(P);
%plot(h0,'k','LineWidth',3);
ylabel('Output')
legend('y(t)')
grid on

figure
subplot(2,1,1)
h1 = plot(t,y,'k','LineWidth',3);
ylabel('Output')
legend('y(t)')
grid on
title("Outpupt and Control Signal for N = " + N + " and Nu = " + Nu)
subplot(2,1,2)
h2 = stairs(t,u,'b','LineWidth',3);
legend('u(t)')
ylabel('Control Signal')
xlabel('Discrete time')
grid on





























