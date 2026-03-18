% Projeto de um filtro de Kalman - Modelo 1 (Linear) - Múltiplas Execuções
% Por: Saulo José Almeida Silva
% Instituição: Universidade Federal de Campina Grande
% Matéria: Sistemas Inteligentes.
% Data: 15/03/2026
% Atualização: 18/03/2026
clc
clear all
pkg load control % Carrega o pacote de controle no Octave

% ---- 1. Sistema Real ----
omega = 1;             % frequência da função seno
dt = 0.01;             % intervalo de amostragem
t = 0:dt:10;           % 10 segundos de dados
N_execucoes = 10;      % Quantidade de simulações para teste de robustez

% Definindo o Modelo Contínuo e Discreto
A = [0 1; -omega^2 0];
B = [0; 0];
grau = length(A);
I = eye(grau);
var_lidas = [1];
H = I(var_lidas, :);
Ad = expm(A*dt);
Hd = H;
Bd = B;

% --- 2. Parâmetros das Covariâncias ---
% (A dedução analítica exata se mantém)
var_w1 = 0.0001;
var_w2 = 0.01;
Q11 = var_w1 * dt + var_w2 * (dt^3)/3;
Q12 = (var_w2 - var_w1 * omega^2) * (dt^2)/2;
Q21 = Q12;
Q22 = var_w2 * dt + var_w1 * (omega^4) * (dt^3)/3;
Qd = [ Q11, Q12 ;
       Q21, Q22 ];
R = 0.01;

% Matriz para armazenar todos os MSEs: linhas = tempo, colunas = execução
todos_mse = zeros(length(t), N_execucoes);

num_sensores = length(var_lidas);

% --- 3. Laço de Múltiplas Execuções (Monte Carlo) ---
for exec = 1:N_execucoes

  % Resetando as condições iniciais para a nova execução
  P = eye(grau);
  x_real = [1; 0];
  x_est  = [1; 0];

  res_real = zeros(length(t), grau);
  res_est  = zeros(length(t), grau);
  z_medido = zeros(length(t), num_sensores);

  % Laço de Simulação e Filtragem para a execução atual
  for i = 1:length(t)
    res_real(i, :) = x_real';

    % Gera medição com novo ruído aleatório
    ruido_sensor = sqrt(R) * randn(num_sensores, 1);
    z = Hd * x_real + ruido_sensor;
    z_medido(i, :) = z';

    % Predição
    x_est_prev = Ad * x_est;
    P_prev     = Ad * P * Ad' + Qd;

    % Correção
    K     = P_prev * Hd' * inv(Hd * P_prev * Hd' + R);
    x_est = x_est_prev + K * (z - Hd * x_est_prev);
    P     = (eye(grau) - K * Hd) * P_prev;

    res_est(i, :) = x_est';

    % Evolui o sistema
    x_real = Ad * x_real;
  endfor

  % Cálculo do MSE da Posição para esta execução específica
  erro_puro = res_real(:,1) - res_est(:,1);
  erro_quadrado = erro_puro.^2;
  mse_progressivo = cumsum(erro_quadrado) ./ (1:length(t))';

  % Guarda o MSE na matriz geral
  todos_mse(:, exec) = mse_progressivo;

endfor

% Calcula a média dos MSEs de todas as execuções
mse_medio_geral = mean(todos_mse, 2);

% --- 4. Plotagem Gráfica ---

figure(1)
% Mostrando os gráficos de rastreio (Apenas da ÚLTIMA execução para referência)
subplot(3,1,1);
plot(t, res_real(:,1), 'b', 'linewidth', 2); hold on;
plot(t, z_medido(:,1), 'r.', 'markersize', 4); hold off;
title(sprintf('1. Posição: Real vs Medição (Apenas Execução %d)', N_execucoes));
ylabel('Amplitude');
legend('Sinal Real', 'Medição Ruidosa', 'location', 'northeast');
grid on;

subplot(3,1,2);
plot(t, res_real(:,1), 'b', 'linewidth', 2); hold on;
plot(t, res_est(:,1), 'g-', 'linewidth', 1.5); hold off;
title('2. Posição: Rastreamento do Filtro');
ylabel('Amplitude');
grid on;

subplot(3,1,3);
plot(t, res_real(:,2), 'b', 'linewidth', 2); hold on;
plot(t, res_est(:,2), 'm-', 'linewidth', 1.5); hold off;
title('3. Velocidade: Rastreamento do Filtro');
xlabel('Tempo (s)');
ylabel('Taxa de Var.');
grid on;

% Gráfico 4: Análise de Robustez - MSE de todas as execuções
figure(2)
plot(t, todos_mse, 'color', [0.7 0.7 0.7], 'linewidth', 1); % Plota as 10 linhas em cinza
hold on;
plot(t, mse_medio_geral, 'k', 'linewidth', 2); % Plota a média em preto grosso
hold off;
title(sprintf('Erro Quadrático Médio (MSE) - %d Execuções', N_execucoes));
xlabel('Tempo (s)');
ylabel('MSE Progressivo');
legend('Execuções Individuais', 'Média Geral', 'location', 'northeast');
grid on;
