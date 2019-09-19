clc;
clear all;

%% leitura de dados e parametrizacoes iniciais
torneira5 = load('torneira5.txt');

% como dito em aula, a entrada ? representada pela 2a coluna do arquivo,
% enquanto y(t) ? dado pela primeira coluna do arquivo.
u_t5 = torneira5(:,2);
y_t5 = torneira5(:,1);
u = -u_t5 + u_t5(1);
y0 = y_t5(1);
Ts = 1;
t = [0:Ts:(length(u_t5)-1)]';

% atraso no tempo theta adquirido por inspecao grafica
theta = 4;
%% plots iniciais sem normaliza
plot(t, u);
hold on;
plot(t, (y_t5-y0)); xlabel('t(s)'); ylabel('Amplitude do sinal');
title('Resposta ao degrau - dados experimentais');
legend('Degrau de entrada u(t)', 'Resposta do sistema y(t)');
grid
hold off
%% normalizando dados
norm_y_t5 = (y_t5 - min(y_t5)) ./ ( max(y_t5) - min(y_t5));
yn =  norm_y_t5;
plot(t, yn); xlabel('t(s)'); ylabel('y(t) normalizado');
title('Resposta ao degrau do sistema normalizada - dados experimentais');
grid
%% Calculo de ganho normalizado e nao normalizado
Knorm = (mean(yn(end-20:end)) - mean(yn(1:20)))/(u(end) - u(1));
K = (mean(y_t5(end-20:end)) - mean(y_t5(1:20)))/(u(end) - u(1));

%% Resposta complementar:


yy = log(abs(1 - yn./(Knorm*u)));
plot(t,yy);
%% calculo do primeiro coeficiente
coef1 = polyfit(t(65:85), yy(65:85), 1);
% ajustes finos
% coef1(1) = coef1(1)-.025;
% coef1(2) = coef1(2)+.45;
% tau1
tau1 = -1/coef1(1);

figure(1); 
subplot(211); plot(t, yy); xlabel('t (s)'); ylabel('w(t)'); axis([60 140 -4 0]);
hold on; plot(t(t>=65), coef1(1)*t(t>=65) + coef1(2), 'm-.', 'LineWidth', 2);
title('Aproximacao da constante de tempo tau1');
subplot(212); plot(t(t>=65), coef1(1)*t(t>=65) + coef1(2), 'm-.', 'LineWidth', 2); 
hold on; plot(t, yy); xlabel('t (s)'); ylabel('w(t)'); axis([65 85 -4 0]);
title('Verificacao da regiao linear para tau1');
grid 
%% calculo do segundo coeficiente
yy2 = log(abs(exp(coef1(2))*exp(-(t)./tau1) - (1 - yn./(Knorm*u))));
plot(t,yy2);
%%
% calculo e ajuste de tau2
coef2 = polyfit(t(65:70),yy2(65:70),1);
% tau2
tau2 = -1/coef2(1);

figure(2); 
subplot(211); plot(t, yy2); xlabel('t (s)'); ylabel('v(t)'); axis([60 140 -10 2]);
hold on; plot(t(t>=65), coef2(1)*t(t>=65) + coef2(2), 'm-.', 'LineWidth', 2); 
title('Aproximacao da constante de tempo tau2');
subplot(212); plot(t(t>=65), coef2(1)*t(t>=65) + coef2(2), 'm-.', 'LineWidth', 2); 
hold on; plot(t, yy2); xlabel('t (s)'); ylabel('v(t)'); axis([65 70 -10 2]);
title('Verificacao da regiao linear para tau2');
grid 

%% resultados
G2a = tf(K, [tau1*tau2  tau1+tau2  1], 'ioDelay', theta);
y2 = lsim(G2a, u, t) + y0;

figure(3); 
plot(t, y2, 'm--','LineWidth',2); xlabel('t (s)'); ylabel('y(t)');
hold on; 
plot(t, y_t5, 'LineWidth',2); xlabel('t (s)'); ylabel('y(t)');

%% dados
K
theta
tau1
tau2
