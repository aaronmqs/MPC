clc, clear, close all;

% Process parameters
tau1 = 20;
tau2 = 50; % it's possible to use negative values for tau
Ts = 1; % Sampling time
z = tf('z',Ts);
s = tf('s');
L = 10;
Ld = ceil(L/Ts);
P = exp(-L*s)/((tau1*s + 1)*(tau2*s + 1));
Pd  = c2d(P,Ts);
[num,den]  = tfdata(Pd,'v');

% CARIMA Model Parameters
a = den;
b = num;
alfa = 0.1; % tuning parameter (robustness)
c = conv([1 -alfa],[1 -alfa]);
delta = [1 -1];
a_til = conv(delta,a);
d = a_til;
[a_til,b] = eqtflength(a_til,b);
[c,d] = eqtflength(c,d); % required for correctly using tf2ss

% Space-State Model
[A,H,B,~] = tf2ss(b,a_til);
A = A';
B = B';
H = H';

[~,~,D,~] = tf2ss(c,d);
D = D';

%% Control Parameters
lambda = 0.001; % Control effort weight in the cost function
delta = 1; % Reference tracking weight in the cost function 
N = 10; % Prediction horizon
Nu = 4; % Control horizon

if Nu > N
    errordlg('Use a Control Horizon not greater than the Prediction Horizon','Invalid Control Horizon.');
    return
end

%% FIR Filter
FFilter = cell(1,Ld+1);
for i = 1:Ld+1
    if i == 1
        FFilter{i} = zeros(size(A*B));
    else
        FFilter{i} = A^(i-2)*B; % *z^-(i-1)
    end
end

% Defining the FIR filter as the coefficients of a discrete time polynomial
FFilter = cell2mat(FFilter);

% Defining the FIR Filter as a discrete time polynomial
% FFilter = (eye(size(A)) - (z^-(Ld))*A^(Ld))*(z*eye(size(A)) - A)\B;

%% Matrix definitions
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

F = cell2mat(F);
G = cell2mat(G);
Qlambda = lambda*eye(size(G,2)); % Control effort weight in the cost function
Qdelta = delta*eye(size(G,1)); % Reference tracking weight in the cost function 

Aux = (G'*Qdelta*G + Qlambda)\G'*Qdelta;
K = Aux(size(B,2),:);
KF = K*F;
KFFIR = KF*FFilter;
Kr = sum(K);

%% Simulation

Tsim = 300;
sim1 = sim('gpc_siso_delay_1_sim');
y = sim1.y;
u = sim1.u;
ref = sim1.ref;
t = 0:Ts:Tsim;

subplot(3,1,1)
y_ol = step(P,t);
plot(t,y_ol,'k',t,ref,'b','LineWidth',3);
ylabel('Output')
title("Open-loop response")
legend('y(t)','Ref')
grid on

subplot(3,1,2)
h1 = plot(t,y,'k',t,ref,'b','LineWidth',3);
ylim([0 1.2])
ylabel('Output')
legend('y(t)','Ref')
grid on
title("Outpupt for N = " + N + " and Nu = " + Nu)
subplot(3,1,3)
h2 = stairs(t,u,'b','LineWidth',3);
title("Control Signal for N = " + N + " and Nu = " + Nu)
legend('u(t)')
ylabel('Control Signal')
xlabel('Discrete time')
grid on