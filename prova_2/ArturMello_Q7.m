clear;
clc
close all;

% G(s) = K/(tau s + 1)

%% numero 7 - leitura de dados

data = load('data_prova_ceai_7.dat');
N = length(data);
t = data(:,1);
u_prbs = data(:,2);
y_out = data(:,3);

Ts = mean(diff(t)); % Intervalo de amostragem dos dados de entrada
% Normalizando saida
% y_out = (y_out - min(y_out))/(max(y_out)-min(y_out));

subplot(2,1,1); plot(t,u_prbs); title('Entrada PRBS do sistema');
xlabel('t(s)'); ylabel('Amplitude');
subplot(2,1,2); plot(t,y_out); title('Resposta do Sistema a entrada PRBS');
xlabel('t(s)'); ylabel('Amplitude');
% Os dados aparentam ser pouco ruidosos

%% definicao do delta de decimacao (secao 12.2.4)

% y_out = sistema superamostrado
% ryy = autocorrelacao de y_out
figure;
subplot(2,1,1);
[tfac1, ryy, lfac1, bfac1] = myccf2(y_out, length(t), 0, 1, 'k');
title('Autocorrelacao r_{y*}'); ylabel('r');
% taumin1 = 85, na funcao acima
taumin1 = tfac1(ryy == min(ryy));

% ryy2 = autocorrelacao de y_out^2
subplot(2,1,2);
[tfac2, ryy2, lfac2, bfac2] = myccf2(y_out.^2, length(t), 0, 1, 'k');
title('Autocorrelacao r_{y*}^2'); ylabel('r');
% taumin2 = 86 na funcao acima
taumin2 = tfac2(ryy2 == min(ryy2));

taumin = min(taumin1, taumin2);
% criterio de selecao de delta para 10 <= taumin/delta <= 20
% definindo que taumin/delta = 14, temos delta:
delta = round(taumin/14); % delta = 6

%% preprocessamento de sinal
% utilizando 70% dos dados para teste, e 30% para validacao
id_size = floor(0.7*N);

t_teste = t(1:id_size);
x_teste = u_prbs(1:id_size);
y_teste = y_out(1:id_size);

t_validacao = t(id_size+1: N);
x_validacao = u_prbs(id_size+1:N);
y_validacao = y_out(id_size+1:N);

delta_k = 1:delta:id_size;

%representacao discreta do sinal
figure; plot(t(delta_k), y_teste(delta_k), 'o'); hold on; plot(t_teste, y_teste);
legend('sinal de teste decimado', 'sinal de teste superamostrado');
xlabel('t(s)'); ylabel('Amplitude');
title('Sinal decimado vs. sinal superamostrado');
% correlacao entrada-saida
figure; myccf2([y_teste(delta_k) x_teste(delta_k)], length(delta_k), 1, 1, 'k');

%% Usando o m?todo de MQ para estimar o modelo do sistema de primeira ordem
% Como visto no livro:
% para sistema de 1a ordem, utilizar os primeiros antecessores 
% para estimar o valor seguinte
% primeiros antecessoers para calculo de parametros
Y = y_teste(1:length(delta_k)-1); 
X = x_teste(1:length(delta_k)-1);
Psi = [Y X]; % matriz de regressores
theta = (Psi' * Psi)\ Psi' * y_teste(2:length(delta_k));

% parametros calculados a partir das informacoes retiradas do livro para um
% sistema de primeira ordem
tau_est = - Ts / (theta(1) - 1)
K_est = (tau_est * theta(2)) / Ts

G = tf([0 K_est],[1 tau_est]);


% validacao dos testes, graficos de dados de validacao e do sistema
% estimado
figure(); 
Y_est = lsim(G, x_validacao, t_validacao);
plot(t_validacao, Y_est);
hold on
plot(t_validacao,y_validacao, 'red')
grid
title('Comparativo entre modelo estimado e dados reais')
legend('modelo estimado', 'resposta do sistema real');

figure(); step(G,t);

% Apesar de o modelo ter apresentado bom formato na estimativa, conformando
% com a resposta lida nos dados do sistema real, ? possivel perceber um
% problema na estimativa. Isso pode ser observado pelo fato de o sistema
% estar inicialmente excitado, ou seja, seu estado nao ? nulo em t=0.