clc;
clear all;

%% leitura de dados e parametrizacoes

planta = load('ENS_25.DAT');
x = planta(:,1);
y = planta(:,2);

% normalizando dados
norm_y = (planta(:,2) - min(planta(:,2))) / ( max(planta(:,2)) - min(planta(:,2)));

plot(x, norm_y, 'ro', 'markersize', 0.5)

% o sistema tem resposta ao degrau semelhante 
% a de um sistema de primeira ordem => H(s) = K/(tau*s + 1)

% pegando a media dos 10 primeiros pontos para definir a amplitude do sinal
%  em y(0):
y_0 = mean(norm_y(1:10));

% pegando a media dos 10 ultimos pontos para definir a amplitude do sinal
%  em y(inf):
y_inf = mean(norm_y((length(norm_y)-9): length(norm_y)));

y_tau = 0.632*(y_inf - y_0) + y_0;
pos_y_tau = find(norm_y >= (y_tau - 0.005) & norm_y <= (y_tau + 0.0001));

% tau sera a media dos valores encontrados onde y_tau eh aproximadamente 63,2% 
% da amplitude maxima do sinal, como calculado acima;

%% construcao do modelo indetificado
tau = x(pos_y_tau(1));
% ajuste de ganho (multiplicado por 0.4 para ajuste fino.
K = (17.05e-3 - 16.34e-3)*1000*0.4;
% modelo = K/((tau*s) + 1)

H_s = tf([0 K],[tau 1]);

% dados de entrada do experimento ENS_25:
% degrau de 16,34 a 17,05 mA ? aplicado.
% opt = stepDataOptions('InputOffset',16.34E-3,'StepAmplitude',17.05E-3);
hold off
plot(x,y, 'ro', 'markersize', 0.5);
hold on
step(H_s, 90:0.01:3700);
legend('Dados experimentais', 'Sistema modelado');
grid

%% validacao do modelo
% invertendo o sinal de validacao para conformar com o modelo proposto
valida = load('ENS_26.DAT');
x = valida(:,1);
y = valida(:,2);
y_inf = mean(y((length(y)-9): length(y)));
% invertendo o sinal para ter validacao do modelo
plot(x,y, 'go', 'markersize', 0.5);
hold on
step(-H_s+y(1)+y_inf, 90:0.01:3700);
legend('Dados experimentais', 'Validacao');
grid

%%
H_s = tf([0 K],[tau 1])
