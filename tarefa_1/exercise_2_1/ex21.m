clc;
clear all;

% parametros iniciais
R = 100;
C = 1E-6;
L = 1E-3;
RC = R*C;
LC = L*C;

s = tf('s'); % funcao de transferencia para realizar rampa unitaria abaixo

% determinando funcao de transferencia para i(s)/v(s)
num_i_v = [0 C 0];
den_i_v = [LC RC 1];
sys_i_v = tf(num_i_v, den_i_v)

polos_i_v = pole(sys_i_v);
% pzmap(sys_i_v) % diagrama de polos e zeros

% determinando funcao de transferencia para vc(s)/v(s)
num_vc_v = [0 0 1/LC];
den_vc_v = [1 (R/L) 1/LC];
sys_vc_v = tf(num_vc_v, den_vc_v)

polos_vc_v = pole(sys_vc_v);
% pzmap(sys_vc_v) % diagrama de polos e zeros


sys = 0;
% option = 0 se sistema em questao e I(s)/V(s)
% option = 1 se sistema em questao e Vc(s)/V(s)
option = 1;
if(option == 0)
    sys = sys_i_v;
elseif(option == 1)
    sys = sys_vc_v;
end

%% exercicio 2.1.3

figure
subplot(311), impulse(sys);   % resposta ao impulso
grid
subplot(312), step(sys);      % resposta ao degrau
grid
subplot(313), step(sys/s);  % resposta a rampa;
grid


%% exercicio 2.1.4

t = 0:0.0001:0.1;
% u = cos(2*pi*60t)
u = 170*cos(377*t);
lsim(sys,u,t);   % u,t define the input signal

%% exercicio 2.1.5
% resposta em frequencia de Vc(s)/V(s)
bode(sys_vc_v);

%% exercicio 2.1.10
t = 0:0.00001:0.001;

u = zeros([length(t) 1]);
spaceState = ss(sys_vc_v);
x0=[1 0];
lsim(spaceState,u,t,x0);   % u,t define the input signal
grid
title('Resposta natural do sistema ao capacitor carregado');

%% exercicio 2.1.11

t = 0:0.00001:0.001;

u = zeros([length(t) 1]);
spaceState = ss(sys_i_v);
x0=[1 0];
lsim(spaceState,u,t,x0);   % u,t define the input signal
grid
title('resposta natural do sistema com corrente no indutor');

