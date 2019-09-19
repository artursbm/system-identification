clc
clear
close all
%% load files
dados = load('dados_tmsd1.txt');

t = dados(:,1);
prbsIn = dados(:,2);
yStar = dados(:,3);
figure;
subplot(2,1,1); plot(t,prbsIn); title('Entrada PRBS do sistema');
xlabel('t(s)'); ylabel('Amplitude');
subplot(2,1,2); plot(t,yStar); title('Resposta do Sistema a entrada PRBS');
xlabel('t(s)'); ylabel('Amplitude');

%% definicao do delta de decimacao (secao 12.2.4)

% yStar = sistema superamostrado
% ryy = autocorrelacao de yStar
figure;
subplot(2,1,1);
[tfac1, ryy, lfac1, bfac1] = myccf2(yStar, length(t), 0, 1, 'k');
title('Autocorrelacao r_{y*}'); ylabel('r');
% taumin1 = 17, na funcao acima
taumin1 = tfac1(ryy == min(ryy));

% ryy2 = autocorrelacao de yStar^2
subplot(2,1,2);
[tfac2, ryy2, lfac2, bfac2] = myccf2(yStar.^2, length(t), 0, 1, 'k');
title('Autocorrelacao r_{y*}^2'); ylabel('r');
% taumin2 = 18 na funcao acima
taumin2 = tfac2(ryy2 == min(ryy2));

taumin = min(taumin1, taumin2);
% criterio de selecao de delta para 10 <= taumin/delta <= 20
% definindo que taumin/delta = 16, temos delta:
delta = round(taumin/16); % delta = 2

%% pre-processamento
N = length(dados);
id_size = floor(0.7*N);

t_teste = t(1:id_size);
x_teste = prbsIn(1:id_size);
y_teste = yStar(1:id_size);

t_validacao = t(id_size+1: N);
x_validacao = prbsIn(id_size+1:N);
y_validacao = yStar(id_size+1:N);

delta_k = 1:delta:id_size;

%representacao discreta do sinal
figure; plot(t(delta_k), y_teste(delta_k), 'o'); hold on; plot(t_teste, y_teste);
legend('sinal de teste decimado', 'sinal de teste superamostrado');
xlabel('t(s)'); ylabel('Amplitude');
title('Sinal decimado vs. sinal superamostrado');
% correlacao entrada-saida
figure; myccf2([y_teste(delta_k) x_teste(delta_k)], length(delta_k), 1, 1, 'k');

%% definindo ordem do sistema usando criterio de Akaike (AIC)
nOrders = 9;
critAic = zeros(1,nOrders);
% abaixo, serao feitas iteracoes para verificar o residuo e obter 
% o valor relacionado ao criterio de akaike, para que seja possivel
% determinar a ordem do modelo estimado que mais o aproxima do sistema real
for i=1:nOrders
    % Primeiro, calculando o residuo:
    yPsi = y_teste(i+1:end);
    [yHat, theta_i] = MinQuadrados(y_teste, x_teste, i);
    residuo = yPsi - yHat;
    % ----- fim de calculo do residuo 
    
    % Calculo do criterio de informacao de Akaike (AIC)
    nAIC = length(yPsi);
    critAic(i) = nAIC * log(var(residuo)) + 2*2*i;
    % ----- fim do calculo de akaike 
end
figure;
plot(critAic)
title('AIC dos dados amostrados')
xlabel('ordem do sistema estimado')
ylabel('AIC')
% conforme estudado, o grafico de AIC denota um joelho em 
% n = 3, ou seja, o modelo estimado tem boa aproximacao de 
% um sistema de 3a ordem

%% Estimativa dos parametros do modelo de 3a ordem
ordem = 3;
[yHat, theta] = MinQuadrados(y_teste, x_teste, ordem);

%% Validacao 1 passo a frente - modelo com a ordem escolhida pelo AIC

% calculo da matriz de regressores para ordem 3 de validacao
PSI = [];
for i = 0:ordem-1
    PSI = [PSI, y_validacao(ordem-i:end-i-1) x_validacao(ordem-i:end-i-1)];
end

yPassoFrente = PSI*theta;

figure()
plot(t_validacao(ordem + 1:end),y_validacao(ordem + 1:end), 'lineWidth', 2)
hold on
plot(t_validacao(ordem + 1:end),yPassoFrente, 'r')
legend;
title('Validacao: Simulacao 1 Passo a Frente');
legend('Sistema Real', 'Modelo Estimado')
xlabel('t')
ylabel('y(t)')
rmsePassoFrente = sqrt(immse(yPassoFrente,y_validacao(ordem+1:end)));

%% Validacao passo livre - modelo com a ordem escolhida pelo AIC

yPassoLivre = zeros(1, length(y_validacao));
yPassoLivre(1:ordem) = y_validacao(1:ordem);
PSI = zeros(length(yPassoLivre), 2*ordem);
PSI(1,:) = [yPassoLivre(1) yPassoLivre(1) x_validacao(1) x_validacao(1)];
PSI(2,:) = [yPassoLivre(1) yPassoLivre(1) x_validacao(1) x_validacao(1)];
PSI(3,:) = [yPassoLivre(2) yPassoLivre(1) x_validacao(2) x_validacao(1)];
for j=ordem+1:length(yPassoLivre)
    PSI(j)=[PSI, yPassoLivre(j-1) yPassoLivre(j-2) x_validacao(j-1) x_validacao(j-2)];
    yPassoLivre(j)=PSI'*theta;
end

figure;
plot(t_validacao, y_validacao, 'lineWidth', 2)
hold on
plot(t_validacao, yPassoLivre);
title('Validacao: Simulacao Livre');
legend('Dados de Validacao', 'Modelo')
