clear;
clc
close all;

% Exercicio 1
%% 1.a) resposta do sistema a PRBS
delay = 10;
Gs = tf([0 2],[20 1],'InputDelay', delay);
Hs = tf([0 2],[20 1]);
b = 9;
N = 500;
Tb = 1;
Ts = 2;

u = prbs(N, b, Tb);
t_u = (0:N-1)*Ts;
y = lsim(Gs, u, t_u);
y0 = lsim(Hs, u, t_u);
figure(1)
% plot com medicao comecando em 0, pois o prbs gerado pela funcao
% acima so gera sinais a partir do tempo t=1 (indice 1 do vetor prbs)
plot(t_u, y)
title('Resposta do sistema Gs com atraso 10 e tempo de amostragem 2 a entrada PRBS')
xlabel('\tau')
ylabel('Y')

%% 1.b) estimacao do atraso pela FCC
[t,r,l,B] = myccf2([y u'],N,1,0,'k');
% plotando a correlacao cruzada da saida pela entrada de acordo com o tempo 
% considerado, como explicado no comentario acima. Este tempo sera 
% iniciado em t=0, para ficar em conformidade.
figure(2)
t_0 = t(1)-1:1:t(length(t))-1;
l_plot = l*ones(length(t_0));
plot(t_0,l_plot,'k:',t_0,-l_plot,'k:',0,1,'k.',0,-1,'k.');
hold on
stem(t_0,r,'k');
hold on
stem(delay/Ts, max(r), 'red'); % pegando valor maximo manualmente para destaque
hold off
xlabel('\tau');
ylabel('r_{uy}');
xlim([-25 25]);
ylim([-1 1]);
title('Funcao de correlacao cruzada da saida pela entrada do sistema Gs');


%% 1.c) MQ para estimacao de parametros
% Metodo dos Minimos Quadrados: para sistema de 1a ordem, pegar os
% primeiros antecessores para estimar o valor seguinte
Y = y0(1:N-1);
X = u(1:N-1)';
Psi = [Y X];
theta = (Psi' * Psi)\ Psi' * y0(2:N);

tau_est = - Ts / (theta(1) - 1)
K_est = (tau_est * theta(2)) / Ts


%% ADICIONANDO RUIDO GAUSSIANO BRANCO BAIXO AO SISTEMA
% 1.d)
power = 0.15;
y_low_noise = y + wgn(N, 1, power);
y0_low_noise = y0 + wgn(N, 1, power);

figure(3)
subplot(2,1,1)
plot(t_u, y_low_noise);
xlabel('\tau')
ylabel('Y LOW NOISE')
title(['Sistema com ruido gaussiano branco de potencia ', num2str(power)])

[t,r,l,B] = myccf2([y_low_noise u'],N,1,0,'k');
subplot(2,1,2)
t_0 = t(1)-1:1:t(length(t))-1;
l_plot = l*ones(length(t_0));
plot(t_0,l_plot,'k:',t_0,-l_plot,'k:',0,1,'k.',0,-1,'k.');
hold on
stem(t_0,r,'k');
hold on
stem(delay/Ts, max(r), 'red'); % pegando valor maximo manualmente para destaque
hold off
xlabel('\tau');
ylabel('r_{uy}');
xlim([-25 25]);
ylim([-1 1]);
title('Funcao de correlacao cruzada da saida pela entrada do sistema');

Y_LOW_NOISE = y0_low_noise(1:N-1);
X = u(1:N-1)';
Psi = [Y_LOW_NOISE X];
theta = (Psi' * Psi)\ Psi' * y0(2:N);

tau_est_low_noise = - Ts / (theta(1) - 1)
K_est_low_noise = (tau_est_low_noise * theta(2)) / Ts

%% ADICIONANDO RUIDO GAUSSIANO BRANCO ALTO AO SISTEMA
power = 5;
y_high_noise=y+wgn(N,1,power);
y0_high_noise = y0 + wgn(N, 1, power);

figure(4)
subplot(2,1,1)
plot(t_u, y_high_noise);
xlabel('\tau')
ylabel('Y HIGH NOISE')
title(['Sistema com ruido gaussiano branco de potencia ', num2str(power)])

subplot(2,1,2)
[t,r,l,B] = myccf2([y_high_noise u'],N,1,0,'k');
t_0 = t(1)-1:1:t(length(t))-1;
l_plot = l*ones(length(t_0));
plot(t_0,l_plot,'k:',t_0,-l_plot,'k:',0,1,'k.',0,-1,'k.');
hold on
stem(t_0,r,'k');
xlabel('\tau');
ylabel('r_{uy}');
xlim([-25 25]);
ylim([-1 1]);
title('Funcao de correlacao cruzada da saida pela entrada do sistema');

Y_HIGH_NOISE = y0_high_noise(1:N-1);
X = u(1:N-1)';
Psi = [Y_HIGH_NOISE X];
theta = (Psi' * Psi)\ Psi' * y0(2:N);

tau_est_high_noise = - Ts / (theta(1) - 1)
K_est_high_noise = (tau_est_high_noise * theta(2)) / Ts

%% Exercicio 2) prbsa02.dat

sys = load('PRBSA02.DAT');
rng = (50:2500);
t = sys(:,1);
u = sys(:,2);
% Normalizando prbs
u = (u - min(u))/(max(u)-min(u));
u = u(rng);
y = sys(:,3);

% filtrando ruido de y
y = medfilt1(y,3);
% Normalizando y
y = (y - min(y))/(max(y)-min(y));
y = y(rng);
N = length(rng);
figure(1)
[t_ccf,r,l,B] = myccf2([y u],N,1,1,'k');
xlabel('\tau');
ylabel('r_{uy}');
title('Funcao de correlacao cruzada da saida pela entrada do sistema');
Ts = mean(diff(t));
Y = y(1:N-1);
X = u(1:N-1);
Psi = [Y X];
theta = (Psi' * Psi)\ Psi' * y(2:N);

tau_est = - Ts / (theta(1) - 1)
K_est = (tau_est * theta(2)) / Ts

Ft = tf([0 K_est],[tau_est 1]);
figure(2)
y_est = lsim(Ft, u, t(rng));
plot(y_est);
hold on;
plot(y, 'red');
grid
title('Comparativo entre dados experimentais e dados estimados');
legend('resposta ao degrau do sistema estimado', 'dados do problema'); 
ylabel('y(\tau)')
xlabel('\tau');
