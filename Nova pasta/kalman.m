% Projeto de um filtro de Kalman
% Por: Saulo José Almeida Silva
% Instituição: Universidade Federal de Campina Grande
% Matéria: Sistemas Inteligentes.
% Data: 15/03/2026
% Atualização: 16/03/2026
clc
clear all
pkg load control % Carrega o pacote de controle no Octave

% ---- 1. Sistema Real
omega = 1;             % frequêcnia da função seno
dt = 0.01;             % intervalo de amostragem
t = 0:dt:10;           % 10 segundos de dados

% DEFININDO O MODELO Modelo contínuo
A = [0 1; -omega^2 0]; % Matriz que gera a senoide
B = [0; 0];            % Não precisa de entrada u, ele oscila sozinho
grau = length(A);

% Definindo a observação do sensor (Matriz H)
I = eye(grau);    % Matriz base de observação
var_lidas = [1];     % Variáveis lidas
H = I(var_lidas, :);

%Modelo discretizado para aplicar em Kalman:
Ad = expm(A*dt);
Hd = H;
Bd = B;

% --- 2. Inicialização das Covariâncias ---
% Necessário conhecer o efeito da discretização no valor de Q
% pois o erro de Q acumula durante o intervalo de amostragem dt
q_var = 0.01; %Valor para aumentar ou diminuir o filtro

%2.1. Discretização da matriz Q (Pelo modelo de ruído da aceleração).
% Calcula como a incerteza se espalha pela posição e velocidade no tempo dt
Qd = q_var * [ (dt^3)/3 , (dt^2)/2 ;
               (dt^2)/2 ,  dt      ];

%2.2. R é escalar (se medimos apenas 1 variável)
R = 0.01;

%2.3. P é a incerteza da nossa estimativa (inicialmente alta = incerto)
P = eye(grau);

% --- 3. Estado Inicial para começar a computagem ---
x_real = [1;0]; %[posição, velocidade]
x_est = [0; 0];


% Resultados do Filtro - Real e estimado.
res_real = [];
res_est = [];
z_medido = [];
erro  = [];

num_sensores  = length(var_lidas);

%Iterações do filtro de kalman
for i = 1:length(t)
  % 1. Salva o estado real ATUAL antes de evoluir
  res_real = [res_real; x_real'];

  % 2. Gera a medição baseada no estado real atual
  ruido_sensor = sqrt(R) * randn(num_sensores, 1);
  z = Hd * x_real + ruido_sensor;
  z_medido = [z_medido; z'];

  % 3. ETAPA KALMAN: Predição (Baseada no x_est anterior)
  x_est_prev = Ad * x_est;
  P_prev     = Ad * P * Ad' + Qd;

  % 4. ETAPA KALMAN: Atualização/Correção
  K     = P_prev * Hd' * inv(Hd * P_prev * Hd' + R);
  x_est = x_est_prev + K * (z - Hd * x_est_prev);
  P     = (eye(grau) - K * Hd) * P_prev;

  % 5. Salva a estimativa corrigida
  res_est = [res_est; x_est'];

  % 6. Evolui o sistema real para o PRÓXIMO passo de tempo
  x_real = Ad * x_real;
endfor

% --- 5. Cálculo do Erro Quadrático Médio (MSE) ---
% Erro instantâneo
erro_puro = res_real(:,1) - res_est(:,1);
% Quadrado do erro
erro_quadrado = erro_puro.^2;
% MSE Progressivo: Média acumulada em cada ponto do tempo
mse_progressivo = cumsum(erro_quadrado) ./ (1:length(t))';

% --- 6. Plotagem ---

% Gráfico 1: Real vs Medição
figure(1)
subplot(2,1,1);
plot(t, res_real(:,1), 'b', 'linewidth', 2); hold on;
plot(t, z_medido(:,1), 'r.', 'markersize', 6); hold off;
title('1. Sistema Real vs Medição Ruidosa');
grid on;

% Gráfico 2: Real vs Estimativa
subplot(2,1,2);
plot(t, res_real(:,1), 'b', 'linewidth', 2); hold on;
plot(t, res_est(:,1), 'g-', 'linewidth', 1.5); hold off;
title('2. Rastreamento do Filtro de Kalman');
legend('Real', 'Kalman');
grid on;

% Gráfico 3: Erro Quadrático Médio (MSE)
figure(2)
plot(t, mse_progressivo, 'm', 'linewidth', 2);
title('3. Erro Quadrático Médio (MSE) Progressivo');
xlabel('Tempo (s)');
ylabel('Valor do MSE');
grid on;
