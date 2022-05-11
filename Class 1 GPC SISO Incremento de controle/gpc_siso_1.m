clc, clear, close all;

% Parâmetros do modelo
% y = (b/a)u + (c/d)e
% CARIMA: a_til*y = b*delta*u + c*e
delta = [1 -1]; % modo da perturbação a ser rejeitada
a = [1 -0.8]; % polos do modelo do processo
a_til = conv(a,delta); 
b = [0 0.4 0.6]; % zeros do modelo do processo
c = 1; % polinômio que indica a robustez

% a = [1 -0.9]; % polos do modelo do processo
% a_til = conv(a,delta); 
% b = [0 0.1]; % zeros do modelo do processo
% c = 1; % polinômio que indica a robustez

d = a_til; % polos do modelo da perturbação
[a_til,b] = eqtflength(a_til,b);
[c,d] = eqtflength(c,d); % necessário para usar o tf2ss

% Parâmetros de controle
N = 3; % Horizonte de predição
Nu = 3; % Horizonte de controle
R = 0.8*eye(N,Nu); % peso do esforço de controle no custo
Q = eye(Nu,N); % peso do erro de seguimento à referência no custo 
r = ones(N,1); % referência

% Space-State Model
[A,H,B,~] = tf2ss(b,a_til);
A = A';
B = B';
H = H';

[~,~,D,~] = tf2ss(c,d);
D = D';

G = cell(N,Nu); 
E = cell(N,1);
F = cell(N,1);

for i = 1:N
        
    for j = 1:Nu
        
        if j == 1
            
            % Encontra dimensão das matrizes para garantir que cada matriz
            % seja um elemento
            dim1 = size(H*A^(i-1)*D);
            if dim1(1)==dim1(2) % garante que não haja elementos repetidos no segundo argumento de num2cell  
                dim1 = dim1(1);
            end
            
            dim2 = size(H*A^i);
            if dim2(1)==dim2(2)  
                dim2 = dim2(1);
            end

            E(i,j) = num2cell(H*A^(i-1)*D,dim1);
            F(i,j) = num2cell(H*A^i,dim2);      

        end
        
        if i >= j
        
            G(i,j) = num2cell(H*A^(abs(i-j))*B);
        
        else
            
            G(i,j) = num2cell(zeros(size(H*B)));
            
        end
        
    end
    
end

E = cell2mat(E);
F = cell2mat(F);
G = cell2mat(G);

M = (G'*Q*G + R)\G'*Q;
K = M(size(B,2),:);
Ke = K*E;
Kf = K*F;
f = Ke + Kf;
Kr = sum(K);

% simulação no simulink
Ts = 1;
Tsim = 20;
sim = sim('gpc_siso_1_sim',Tsim);
y = sim.y;
delta_u = sim.delta_u;
t = 0:Ts:Tsim;

plot(t,y,'LineWidth',3);
hold on
stairs(t,delta_u,'LineWidth',3)
xlabel('Discrete label')
legend('y(t)','\Deltau(t)')
grid on
















