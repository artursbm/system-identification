clc;
clear all;

to = 0;
tff = 300;
ode_o = [1;0.8;0.5];

[t,output] = ode45(@food_chain_ode_sim, [to tff], ode_o);

% plot(t,output)
% grid
% title('Populacoes normalizadas de Predador Generalista, Predador especialista e Presa ao longo do tempo')
% xlabel('tempo')
% ylabel('Populacao normalizada')
% legend('Presa', 'Predador especialista', 'Predador generalista', 'Location', 'North')


plot(output(:,1), output(:,2))
title('Evolucao de Presas em relacao a Predadores Especialistas')
xlabel('Presas')
ylabel('Predador Especialista')

% Dados retirados da simulacao realizada acima.
x_ini=[0.2989 0.260 0.227 0.199 0.179 0.165 0.159 0.158 0.1630];
% 
% 
for i=(9:100)
    out = model_2(x_ini, i);
    x_ini(i) = out;
end

plot(x_ini, '-')
hold on
% 0.25 de ajuste no modelo para se chegar a um modelo satisfatorio
plot((0:100),(output(160:260,1)+0.25))
title('Relacao entre modelos 1 e 2 para intervalo de tempo entre 160 e 260')
xlabel('Tempo t')
ylabel('Populacao normalizada')
legend('Modelo 1', 'Modelo 2', 'Location', 'North')

