% Projeto de um filtro de Kalman - Modelo 4 (EKF Adaptativo: y, y_dot, omega)
% Por: Saulo José Almeida Silva
% Atualização: 18/03/2026
clc; clear all; close all;
pkg load control;

% ---- 1. Configurações da Simulação Real ----
A_real = 1.5;           % Amplitude real
omega_real = 2.0;       % Frequência real (rad/s)
phi0_real = 0;          % Fase inicial
dt = 0.01;              % Período de amostragem
TMAX = 10;              % Tempo máximo
t = (0:dt:TMAX)';       % Vetor de tempo
N_execucoes = 10;

% ---- 2. Parâmetros do Filtro ----
grau = 3;               % Estado: [y; y_dot; omega]
H = [1, 0, 0];          % Observação linear
R = 0.05;               % Covariância da medição

% Incerteza da Frequência (usada para montar Qk)
sig_omega_sq = 0.005;

% ---- 3. Simulação de Monte Carlo ----
mse_y      = zeros(length(t), N_execucoes);
mse_ydot   = zeros(length(t), N_execucoes);
mse_omega  = zeros(length(t), N_execucoes);

% Vetores reais para comparação (Gabarito)
y_real_vetor     = A_real * sin(omega_real * t + phi0_real);
ydot_real_vetor  = A_real * omega_real * cos(omega_real * t + phi0_real);
omega_real_vetor = omega_real * ones(length(t), 1);

for exec = 1:N_execucoes
    % Condições Iniciais: [y; y_dot; omega]
    % Chute inicial: Começa com frequência errada (1.0) para testar convergência
    x_est = [y_real_vetor(1); ydot_real_vetor(1); 1.0];
    P = eye(grau) * 0.5;

    est_param = zeros(length(t), grau);

    for i = 1:length(t)
        % --- GERAÇÃO DA MEDIÇÃO ---
        z = y_real_vetor(i) + sqrt(R) * randn();

        % Estados atuais para as Jacobianas
        y_hat  = x_est(1);
        yd_hat = x_est(2);
        w_hat  = x_est(3);

        % --- EKF: PREDIÇÃO NÃO-LINEAR ---
        x_pred = [ y_hat + yd_hat * dt;
                   yd_hat - (w_hat^2 * y_hat) * dt;
                   w_hat ];

        % --- JACOBIANAS ADAPTATIVAS (Fk e Qk) ---
        Fk = [ 1,             dt,     0;
              -w_hat^2 * dt,  1,     -2 * w_hat * y_hat * dt;
               0,             0,      1 ];

        Qk = sig_omega_sq * [ 0, 0, 0;
            0, (4 * w_hat^2 * y_hat^2 * (dt^3)/3), (-2 * w_hat * y_hat * (dt^2)/2);
            0, (-2 * w_hat * y_hat * (dt^2)/2), dt ];

        P_pred = Fk * P * Fk' + Qk;

        % --- EKF: CORREÇÃO ---
        y_tilde = z - H * x_pred;
        S = H * P_pred * H' + R;
        K = (P_pred * H') / S;

        x_est = x_pred + K * y_tilde;
        P = (eye(grau) - K * H) * P_pred;

        est_param(i, :) = x_est';
    end

    % --- CÁLCULO DOS ERROS (Execução Atual) ---
    err_y     = y_real_vetor - est_param(:, 1);
    err_ydot  = ydot_real_vetor - est_param(:, 2);
    err_omega = omega_real_vetor - est_param(:, 3);

    mse_y(:, exec)     = cumsum(err_y.^2) ./ (1:length(t))';
    mse_ydot(:, exec)  = cumsum(err_ydot.^2) ./ (1:length(t))';
    mse_omega(:, exec) = cumsum(err_omega.^2) ./ (1:length(t))';
end

% Médias de Monte Carlo
mse_medio_y     = mean(mse_y, 2);
mse_medio_ydot  = mean(mse_ydot, 2);
mse_medio_omega = mean(mse_omega, 2);

% ---- 4. Plotagem ----

% Figura 1: Rastreamento e Convergência (Última Execução)
figure(1);
subplot(3,1,1);
plot(t, y_real_vetor, 'b', 'LineWidth', 1.5); hold on;
plot(t, est_param(:,1), 'r--', 'LineWidth', 1);
title('Modelo 4: Rastreamento do Sinal (y)'); ylabel('Amplitude'); legend('Real','EKF'); grid on;

subplot(3,1,2);
plot(t, ydot_real_vetor, 'b', 'LineWidth', 1.5); hold on;
plot(t, est_param(:,2), 'g--', 'LineWidth', 1);
title('Modelo 4: Derivada (\dot{y})'); ylabel('Taxa'); grid on;

subplot(3,1,3);
plot(t, omega_real_vetor, 'k--', 'LineWidth', 1.5); hold on;
plot(t, est_param(:,3), 'm', 'LineWidth', 1.5);
title('Modelo 4: Convergência da Frequência (\omega)'); ylabel('rad/s'); grid on;

% Figura 2: Análise de Monte Carlo (MSE do Sinal y)
figure(2);
plot(t, mse_y, 'Color', [0.7 0.7 0.7], 'LineWidth', 1); hold on;
plot(t, mse_medio_y, 'k', 'LineWidth', 2);
title(sprintf('Modelo 4: MSE do Sinal - %d Execuções', N_execucoes));
xlabel('Tempo (s)'); ylabel('MSE Progressivo');
legend('Execuções Individuais', 'Média Geral'); grid on;

% Figura 3: MSE Médio dos outros estados
figure(3);
subplot(2,1,1);
plot(t, mse_medio_ydot, 'g', 'LineWidth', 2);
title('MSE Médio: Derivada (\dot{y})'); ylabel('MSE'); grid on;

subplot(2,1,2);
plot(t, mse_medio_omega, 'r', 'LineWidth', 2);
title('MSE Médio: Frequência (\omega)'); ylabel('MSE'); xlabel('Tempo (s)'); grid on;d
