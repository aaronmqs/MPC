clc, clear, close all;

% Diferente do programa da aula 1, neste o sinal de controle é calculado
% diretamente, em vez do esforço de controle
% Parâmetros do modelo
% y = (b/a)u + (c/d)e
% CARIMA: a_til*y = b*delta*u + c*e
delta = [1 -1]; % modo da perturbação a ser rejeitada
a = [1 -0.8]; % polos do modelo do processo
a_til = conv(a,delta); 
b = [0 0.4 0.6]; % zeros do modelo do processo
alfa = 0.1;
c = conv([1 -alfa],[1 -alfa]); % polinômio que indica a robustez
% c = 1;
d = a_til; % polos do modelo da perturbação
[a,b] = eqtflength(a,b);
[c,d] = eqtflength(c,d); % necessário para usar o tf2ss

% Space-State Model
[A,H,D,~] = tf2ss(c,d);
A = A';
H = H';
D = D';

[Ab,Hb,Bb,~] = tf2ss(b,a);
Ab = Ab';
Bb = Bb';
Hb = Hb';

N = 2; % prediction horizon
Nu = 1; % Horizonte de controle

if Nu > N
    errordlg('Use a Control Horizon smaller than the Prediction Horizon','Invalid Control Horizon.');
    return
end

G = cell(N); 
E = cell(N,1);
F = cell(N,1);
Fb = F;

for i = 1:N
    F{i,1} = H*A^i;
    E{i,1} = H*A^(i-1)*D;
    Fb{i,1} = Hb*Ab^i;
    
    for j = 1:N
        if i < j
            G{i,j} = zeros(size(Hb*Bb,1));
        else
            G{i,j} = Hb*Ab^(i-j)*Bb;
        end
    end

end

% Y = G*U + Fb*x(t) + F*p(t) + E*e(t), com U = [u(t) u(t+1) ... u(t + N - 1)]
% Para o caso em que Nu < N, deve-se observar que as (N - Nu) últimas
% entradas do vetor U são todas idênticas a u(t + Nu - 1) e, portanto,
% deve-se representar a matriz G com somente Nu colunas, sendo a última
% delas igual à soma das (N - Nu + 1) últimas colunas da matriz original G.

% Prepares the Columnn Nu with the sum of the other columns 
for i = 1:N
    for j = Nu+1:N
        G{i,Nu} = G{i,Nu} + G{i,j};
    end
end

Gaux = cell(N,Nu);

% Builds an auxiliar matrix with only Nu columns from G
for i = 1:N
    for j = 1:Nu
        Gaux{i,j} = G{i,j};
    end
end
G = Gaux;

E = cell2mat(E);
F = cell2mat(F);
Fb = cell2mat(Fb);
G = cell2mat(G);

% Transformation from Delta_U to U
Bo = cell(Nu,1);
M = cell(Nu);

for i = 1:Nu

    if i == 1

        Bo{i,1} = eye(size(Bb,2));

    else

        Bo{i,1} = zeros(size(Bb,2));

    end

    for j = 1:Nu

        if i == j

            M{i,j} = eye(size(Bb,2));

        elseif i == j+1

            M{i,j} = -eye(size(Bb,2));

        else

            M{i,j} = zeros(size(Bb,2));

        end
    end

end

M = cell2mat(M);
Bo = cell2mat(Bo);

% Parâmetros de controle
lambda = 0.8;
delta = 1;
Qlambda = lambda*eye(size(G,2)); % peso do esforço de controle no custo (Q_lambda)
Qdelta = delta*eye(size(G,1)); % peso do erro de seguimento à referência no custo 

Aux = (G'*Qdelta*G + M'*Qlambda*M)\G'*Qdelta;
K = Aux(size(Bb,2),:);
Ke = K*E;
Kf = K*F;
Kfb = K*Fb;
f = Ke + Kf + Kfb;
Kr = sum(K);

% This gain appears because of the transformation from Delta_U to U
Aux0 = (G'*Qdelta*G + M'*Qlambda*M)\M'*Qlambda*Bo;
Ko = Aux0(size(Bb,2),:);

% simulação no simulink
Ts = 1;

% Bloco V
z = tf('z',Ts);
V = K*(F - E*H)/(z*eye(size(A,1)) - A + D*H)*D + K*E;

Tsim = 30;
sim1 = sim('gpc_siso_2_sim');
y = sim1.y;
u = sim1.u;
y1 = sim1.y1;
u1 = sim1.u1;
t = 0:Ts:Tsim;

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













