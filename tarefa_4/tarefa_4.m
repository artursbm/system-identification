clear all
clc

%% Exercicio 4.8

j=1000;
% gerando ruido branco com distribuicao gaussiana e variancia 1;
e = wgn(j);

figure(1)
plot(e)

u=zeros(1,j);

for k=4:j
    u(k) = 0.9*e(k-1) + 0.8*e(k-2) + 0.7*e(k-3) + e(k);
end

figure(2)
[t,r,l,B] = myccf2(u', 5, 0, 1, 'k');
title('Autocorrelacao do sinal u(k)')
grid

%% Exercicio 4.15

N = 189;
b = 6;
m = 1;

y = prbs(N,b,m);
subplot(2,1,1);
stairs(y);
title('PRBS de sequencia m com 6 bits de resolucao e 3 `periodos` de amostragem'); 

subplot(2,1,2);
res = myccf2(y', N, 1, 1, 'k');
title('FAC do sinal PRBS de sequencia m de 6 bits e 3 `periodos` (189 pontos)');
grid

%% Exercicio 4.16

Tb = [1 2 5 10];
b = 6;
N = (2^b - 1)*3*Tb; % 189*Tb
m=Tb;
j=0;
for i=1:length(Tb)
    j=j+1;
    y = prbs(N(i),b,m(i));
    
    figure(j)
    subplot(2,1,1);
    stairs(y);
    title(['PRBS de sequencia m com 6 bits de resolucao e Tb =' num2str(Tb(i))]); 
    subplot(2,1,2);
    [t,r,l,B]=myccf2(y', N(i), 1, 1, 'k');
    xlim([-200 200])
    title(['FAC do sinal PRBS de sequencia m de 6 bits e Tb = ' num2str(Tb(i))]);
    grid
    
    
end

%% Exercicio 4.20

H = tf([0 1],[1000 1]);
Tb = [1 100 1000 10000];
b = 12;
% multiplica-se N por 4 para que o PRBS seja amostrado no periodo
% suficiente para que o sistema consiga ter alguma resposta relevante
% quando Tb >= 10000. 
N = (2^b - 1)*3;
m=Tb;
for t=1:length(Tb)
    
    u = prbs(N, b, m(t));
    figure(t)
    lsim(H, u, 1:1:N);
    xlim([0 1E4])
    title(['Resposta do sistema ao PRBS com Tb=' num2str(Tb(t))]);

end
