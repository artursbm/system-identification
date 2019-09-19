clc;
clear all;

data = load('DS100g_prato_1.txt');
data_x = data(:,1);
% norm_data_y = data(:,2);
% uns = ones(length(data(:,2)));

% normalizando a entrada de dados
norm_data_y = (data(:,2) - min(data(:,2))) / ( max(data(:,2)) - min(data(:,2)));
uns = (ones(length(data(:,2))) - min(data(:,2))) / ( max(data(:,2)) - min(data(:,2)));

% transformando as amostras para o dominio de tempo a partir da amostragem
% de Ts = 5ms;
samples = data_x;         
Fs = 1/0.005;             
t = samples/Fs; 

t = t(t>=0);
norm_y = norm_data_y(250:length(norm_data_y));
% segunda ordem subamortecido com coeficiente baixo
plot(t, norm_y, 'black');
% findpeaks(norm_y)
yPeaks = findpeaks(norm_y);
% hold on;
% plot(t, uns);
% yaxis = 0:0.1:1.0;
% plot(ones(length(yaxis))*0, yaxis);
% plot(ones(length(yaxis))*0.014, yaxis);

grid

% ganho: K por ajuste fino
% frequencia natural: wn
% numero de ciclos: N
% coeficiente de amortecimento: zeta
K=0.35;
N = length(yPeaks);
T = 3E-4;
wn = (2*pi)/T;
zeta = 0.6/N;

%% resultados
fprintf('frequencia natural: %.2f rad/s, coeficiente de amortecimento: %.5f;\n',wn, zeta);
hold on
Hs = tf([0 0 (K*wn^2)],[1 (2*zeta*wn) wn^2])
step(Hs+norm_y(1), 0.02)
title('Step response')
legend('Dados experimentais', 'Modelo estimado');


%% validacao
data = load('DS100g_prato_2.txt');
data_x = data(:,1);
% norm_data_y = data(:,2);
% uns = ones(length(data(:,2)));

% normalizando a entrada de dados
norm_data_y = (data(:,2) - min(data(:,2))) / ( max(data(:,2)) - min(data(:,2)));
uns = (ones(length(data(:,2))) - min(data(:,2))) / ( max(data(:,2)) - min(data(:,2)));

samples = data_x;         
Fs = 1/0.005;             
t = samples/Fs; 

t = t(t>=0);
norm_y = norm_data_y(250:length(norm_data_y));
% segunda ordem subamortecido com coeficiente baixo
plot(t, norm_y, 'black');
hold on
Hs = tf([0 0 (K*wn^2)],[1 (2*zeta*wn) wn^2]);
step(Hs+norm_y(1), 0.02)
title('Step response - Validacao do sistema')
legend('Dados experimentais', 'Modelo 2 estimado');
