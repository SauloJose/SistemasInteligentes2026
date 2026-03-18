% Projeto de um filtro de Kalman - Modelo 2 (EKF) - Monte Carlo
% Por: Saulo José Almeida Silva
% Atualização: 18/03/2026
clc; clear all; close all;
pkg load control;

% ---- 1. Configurações da Simulação Real ----
A_real = 1.5;           % Amplitude real
omega_real = 1.2;       % Frequência real (rad/s)
phi0_real = 0;          % Fase inicial
dt = 0.01;              % Período de amostragem
TMAX = 10;              % Tempo máximo da amostragem
t = (0:dt:TMAX)';       % Vetor de tempo
N_execucoes = 10;

% ---- 2. Parâmetros do Filtro (Modelo 2) ----
% Estado: [Fase; Frequência; Amplitude]
grau = 3;

% Matriz de Transição Discreta (Exata pois F^2 = 0)
F = [1 dt 0;
     0 1  0;
     0 0  1];

% Covariância do Ruído de Processo (Qc[1,1] = 0 conforme pedido)
sig_omega = 0.001; % Incerteza na frequência
sig_A = 0.02;     % Incerteza na amplitude

Qk = [ (sig_omega^2)*(dt^3)/3, (sig_omega^2)*(dt^2)/2, 0;
       (sig_omega^2)*(dt^2)/2, (sig_omega^2)*dt,       0;
       0,                      0,               (sig_A^2)*dt ];

% Covariância da Medição (Escalar)
R = 0.05;

% ---- 3. Simulação de Monte Carlo ----
todos_mse = zeros(length(t), N_execucoes);

% Matrizes para guardar o MSE de cada estado
mse_fase = zeros(length(t), N_execucoes);
mse_freq = zeros(length(t), N_execucoes);
mse_amp  = zeros(length(t), N_execucoes);

% Vetores reais dos estados para comparação ao longo do tempo
fase_real_vetor = omega_real * t + phi0_real; % A fase real varia no tempo
freq_real_vetor = omega_real * ones(length(t), 1); % Constante
amp_real_vetor  = A_real * ones(length(t), 1);     % Constante

for exec = 1:N_execucoes
    % Condições Iniciais do Filtro (Estimativa inicial com erro)
    x_est = [0.1; 1.0; 1.0]; % Começa com valores diferentes do real [fase, fre. ang, Amplitude]
    P = eye(grau) * 0.1;

    res_real = A_real * sin(omega_real * t + phi0_real);
    res_est  = zeros(length(t), 1);
    est_param = zeros(length(t), grau); % Para salvar phi, omega, A

    for i = 1:length(t)
        % --- GERAÇÃO DA MEDIÇÃO REAL ---
        z = res_real(i) + sqrt(R) * randn();

        % --- EKF: ETAPA DE PREDIÇÃO ---
        % Transição é linear: x = F*x
        x_pred = F * x_est;
        P_pred = F * P * F' + Qk;

        % --- EKF: ETAPA DE CORREÇÃO ---
        % 1. Jacobiana da Medição H_k = [A*cos(phi), 0, sin(phi)]
        phi_p = x_pred(1);
        w_p   = x_pred(2);
        A_p   = x_pred(3);

        Hk = [A_p * cos(phi_p), 0, sin(phi_p)];

        % 2. Inovação (Usa a função não-linear h(x))
        h_x = A_p * sin(phi_p);
        y_tilde = z - h_x;

        % 3. Ganho de Kalman
        S = Hk * P_pred * Hk' + R;
        K = (P_pred * Hk') / S;

        % 4. Atualização
        x_est = x_pred + K * y_tilde;
        P = (eye(grau) - K * Hk) * P_pred;

        % Salva resultados
        res_est(i) = x_est(3) * sin(x_est(1)); % Sinal reconstruído
        est_param(i, :) = x_est';
    end

    % Cálculo do MSE do sinal (como já estava)
    erro_puro = res_real - res_est;
    todos_mse(:, exec) = cumsum(erro_puro.^2) ./ (1:length(t))';

    % Cálculo do MSE do Individual (como já estava)
    erro_fase = fase_real_vetor - est_param(:, 1);
    erro_freq = freq_real_vetor - est_param(:, 2);
    erro_amp  = amp_real_vetor  - est_param(:, 3);

    mse_fase(:, exec) = cumsum(erro_fase.^2) ./ (1:length(t))';
    mse_freq(:, exec) = cumsum(erro_freq.^2) ./ (1:length(t))';
    mse_amp(:, exec)  = cumsum(erro_amp.^2)  ./ (1:length(t))';
end

% Calcula a média do MSE ao longo das execuções de Monte Carlo
mse_medio_geral = mean(todos_mse, 2);
mse_medio_fase  = mean(mse_fase, 2);
mse_medio_freq  = mean(mse_freq, 2);
mse_medio_amp   = mean(mse_amp, 2);

% ---- 4. Plotagem ----
figure(1);
subplot(3,1,1);
plot(t, res_real, 'b', 'LineWidth', 1.5); hold on;
plot(t, res_est, 'r--', 'LineWidth', 1);
title('Rastreamento do Sinal (Última Execução)');
legend('Real', 'EKF'); grid on;

subplot(3,1,2);
plot(t, est_param(:,2), 'm', 'LineWidth', 1.5); hold on;
line([t(1) t(end)], [omega_real omega_real], 'Color', 'k', 'LineStyle', '--');
title('Convergência da Frequência (\omega)');
legend('Estimado', 'Real'); grid on;

subplot(3,1,3);
plot(t, est_param(:,3), 'm', 'LineWidth', 1.5); hold on;
line([t(1) t(end)], [A_real A_real], 'Color', 'k', 'LineStyle', '--');
title('Convergência da Amplitude (A)');
legend('Estimado', 'Real'); grid on;

figure(2);
plot(t, todos_mse, 'Color', [0.8 0.8 0.8]); hold on;
plot(t, mse_medio_geral, 'k', 'LineWidth', 2);
title('Análise de Monte Carlo: MSE do Sinal');
xlabel('Tempo (s)'); ylabel('MSE');
grid on;

% ---- Gráficos do MSE dos Estados Individuais ----
figure(3);
subplot(3,1,1);
plot(t, mse_medio_fase, 'b', 'LineWidth', 2);
title('Análise de Monte Carlo: MSE Médio da Fase');
xlabel('Tempo (s)'); ylabel('MSE');
grid on;

subplot(3,1,2);
plot(t, mse_medio_freq, 'r', 'LineWidth', 2);
title('Análise de Monte Carlo: MSE Médio da Frequência');
xlabel('Tempo (s)'); ylabel('MSE');
grid on;

subplot(3,1,3);
plot(t, mse_medio_amp, 'g', 'LineWidth', 2);
title('Análise de Monte Carlo: MSE Médio da Amplitude');
xlabel('Tempo (s)'); ylabel('MSE');
grid on;
