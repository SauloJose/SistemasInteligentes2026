% Projeto de um filtro de Kalman
% Por: Saulo José Almeida Silva
% Instituição: Universidade Federal de Campina Grande
% Matéria: Sistemas Inteligentes.
% Data: 16/03/2026
% Atualização: 17/03/2026

clc; clear all; close all;
pkg load control % Carrega o pacote de controle no Octave

% ---- 1. Parâmetros Gerais
omega_real = 1;        % Frequência real (desconhecida para o filtro)
dt = 0.01;             % Intervalo de amostragem
t = 0:dt:10;           % 10 segundos de simulação

% --- 2. Sistema Real (Planta) ---
% Continuamos simulando o mundo real de forma linear exata para gerar os dados
A_real = [0 1; -omega_real^2 0];
Ad_real = expm(A_real*dt);
H_real = [1 0];

% --- 3. Configuração do EKF ---
% O novo vetor de estados tem 3 variáveis: [posição; velocidade; frequencia_quadrada]
grau_ekf = 3;
Hd = [1 0 0]; % O sensor continua lendo apenas a posição (x1)

% Inicialização das Covariâncias
q_var = 0.01;
Qd = zeros(grau_ekf, grau_ekf);

% Ruído para posição e velocidade
Qd(1:2, 1:2) = q_var * [ (dt^3)/3 , (dt^2)/2 ;
                         (dt^2)/2 ,  dt      ];
% Ruído para a frequência quadrada (permite que o filtro ajuste esse parâmetro)
Qd(3,3) = 1e-4;

R = 0.01;
P = eye(grau_ekf); % Incerteza inicial

% --- 4. Estados Iniciais ---
x_real = [1; 0];          % Real: [posição, velocidade]
% O Filtro começa achando que a frequência quadrada é 0.5 (errado, o real é 1.0)
x_est = [0; 0; 0.5];      % Est: [posição, velocidade, mu]

% Armazenamento
res_real = []; res_est = []; z_medido = [];

% --- 5. Laço do EKF ---
for i = 1:length(t)
  % 1. Salva o estado real
  res_real = [res_real; x_real'];

  % 2. Medição Ruidosa
  ruido_sensor = sqrt(R) * randn(1, 1);
  z = H_real * x_real + ruido_sensor;
  z_medido = [z_medido; z'];

  % 3. ETAPA KALMAN: Predição Não-Linear (Euler Simples Vetorizado)

  % Vetor de derivadas f(x) = [velocidade; aceleração; derivada_do_parametro]
  f_x = [ x_est(2);
         -x_est(3) * x_est(1);
          0 ];

  % Método de Euler Simples em uma única linha matricial
  x_est_prev = x_est + f_x * dt;

  % Jacobiano do sistema contínuo avaliado no estado atual
  F_continuo = [ 0,         1,  0;
                -x_est(3),  0, -x_est(1);
                 0,         0,  0 ];

  % Jacobiano discretizado via Exponencial de Matriz (Mais preciso)
  Fd = expm(F_continuo * dt); %e^(F_continuo*dt)

  % Propaga a covariância
  P_prev = Fd * P * Fd' + Qd;

  % 4. ETAPA KALMAN: Atualização/Correção
  Inovacao = z - Hd * x_est_prev;
  S = Hd * P_prev * Hd' + R;
  K = P_prev * Hd' / S; % Usando '/' no lugar de 'inv()' por boas práticas numéricas

  x_est = x_est_prev + K * Inovacao;
  P = (eye(grau_ekf) - K * Hd) * P_prev;

  % 5. Salva a estimativa corrigida
  res_est = [res_est; x_est'];

  % 6. Evolui o sistema real
  x_real = Ad_real * x_real;
endfor

% --- 5. Cálculo do Erro Quadrático Médio (MSE) ---
% Erro instantâneo
erro_puro = res_real(:,1) - res_est(:,1);
% Quadrado do erro
erro_quadrado = erro_puro.^2;
% MSE Progressivo: Média acumulada em cada ponto do tempo
mse_progressivo = cumsum(erro_quadrado) ./ (1:length(t))';


% --- 6. Plotagem ---
figure('Name', 'Resultados EKF');

figure(1)
subplot(2,1,1);
plot(t, res_real(:,1), 'b', 'linewidth', 2); hold on;
plot(t, z_medido(:,1), 'r.', 'markersize', 6); hold off;
title('1. Sistema Real vs Medição Ruidosa (Posição)'); grid on;

subplot(2,1,2);
plot(t, res_real(:,1), 'b', 'linewidth', 2); hold on;
plot(t, res_est(:,1), 'g-', 'linewidth', 1.5); hold off;
title('2. Rastreamento da Posição (EKF)'); legend('Real', 'EKF'); grid on;

figure(2)
% Linha preta tracejada no valor 1.0 (que é 1^2)
plot(t, ones(size(t))*omega_real^2, 'k--', 'linewidth', 2); hold on;
plot(t, res_est(:,3), 'm-', 'linewidth', 2); hold off;
title('3. Estimação do Parâmetro Desconhecido (\mu = \omega^2)');
legend('Valor Real', 'Estimativa EKF'); xlabel('Tempo (s)'); grid on;

% Gráfico 3: Erro Quadrático Médio (MSE)
figure(3)
plot(t, mse_progressivo, 'm', 'linewidth', 2);
title('3. Erro Quadrático Médio (MSE) Progressivo');
xlabel('Tempo (s)');
ylabel('Valor do MSE');
grid on;
