clear
clc
close all

% SISTEMA DE BOMBEAMENTO DE AGUA

%% Espaco de estados
rho = 1;
A1 = 20;
A2 = 15;
g = 9.8;
L = 10;
A3 = 2;
% x1 = H1, x2 = H2, x3 = W3
A = [0 0 -1/(rho*A1); 0 0 1/(rho*A2); (g*rho*A3)/L, -(g*rho*A3)/L, 0];
B = [0; -1/(rho*A2); 0];

C = [1 0 0; 0 1 0];
D = [0; 0];

SYS = ss(A,B,C,D);
SYS_inv = ss(-A,-B,C,D);

t = 0:0.01:50;

u_step = t<=t(floor(length(t)/2));
% lsim(SYS, u, t);
figure(); lsim(SYS_inv, u_step, t);
% figure(); step(SYS,t);

youts = lsim(SYS_inv, u_step, t);

%% estimando FT

[num, den] = ss2tf(-A,-B,C,D);

H1 = tf(num(1,:), den)
figure(); bode(H1);

H2 = tf(num(2,:), den)
figure(); bode(H2);

%% inserindo ruido e calculando sistema de primeira ordem
y_h1 = youts(:,1);
y_h2 = youts(:,2);

plot(t,y_h1);
title('Resposta ao degrau do sistema estimado para o tanque 1');
xlabel('tempo(s)');
ylabel('Amplitude');
% visivelmente o atraso de tempo em H1 ? 1.35

N = length(youts);
b = floor(log2(N)+1);
u_prbs = prbs(N,b,10);

% gerando sinal de H1 com ruido branco de 20dB
y_noise = awgn(y_h1,20);

figure(); plot(t,y_noise);
title('Resposta ao degrau do sistema do tanque 1 com ruido branco de 20dB');
xlabel('tempo(s)');
ylabel('Amplitude');

% estimando ganho pela media no periodo de estado permanente (t>=25)
K = mean(y_noise(t>=20));
y_tau = 0.632*K;
all_tau = t(y_noise>=y_tau-0.0005 & y_noise <=y_tau+0.0005);
tau = mean(all_tau);

num = [0 K];
den = [tau 1];

H_est = tf(num, den);

figure();
step(H_est, t, 'red')
hold on
plot(t,y_h1);
title('Resposta ao degrau do modelo real vs. modelo estimado');
xlabel('tempo(s)');
ylabel('Amplitude');

%% Sistema integrador

theta = 0.54;
delta_u = 1;
delta_t_u = 25;

delta_y = mean(y_noise(t>=25))
K = delta_y/(delta_t_u * delta_u)

t_est = t;
f_est = ones(1,length(t_est));
f_est(1:301) = 0;
f_est(2501:length(t_est)) = delta_u*K*delta_t_u;
for i=302:2500
    f_est(i) = K*t_est(i)*delta_u;
end

G_est = tf([0 K],[1 0], 'InputDelay', theta);

figure(); 
plot(t_est, f_est, 'g');
hold on
plot(t,y_h1, 'r');
hold on
lsim(G_est, u_step, t);
title('Resposta ao pulso do sistema real vs. modelo estimado');
xlabel('tempo(s)');
ylabel('Amplitude');
legend('Modelo estimado','Sistema real', 'Sistema integrador');
