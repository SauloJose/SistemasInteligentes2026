% Projeto de um filtro de Kalman - Modelo 5 (EKF: y, y_dot, z=omega^2) - Monte Carlo
% Por: Saulo José Almeida Silva
% Atualização: 18/03/2026
clc; clear all; close all;
pkg load control;

% ---- 1. Configurações da Simulação Real ----
A_real = 1.5;            % Amplitude real
omega_real = 1.2;       % Frequência real (rad/s)
z_real = omega_real^2;  % Frequência ao quadrado real
dt = 0.01;              % Período de amostragem
TMAX = 10;              % Tempo máximo
t = (0:dt:TMAX)';       % Vetor de tempo
N_execucoes = 10;

% ---- 2. Parâmetros do Filtro ----
grau = 3;               % Estado: [y; y_dot; z]
H = [1, 0, 0];          % Observação linear
R = 0.1;                % Covariância da medição (Ruído do sensor)
sig_z_sq = 0.001;         % Incerteza do processo sobre z

% ---- 3. Simulação de Monte Carlo ----
todos_mse = zeros(length(t), N_execucoes);
mse_y     = zeros(length(t), N_execucoes);
mse_ydot  = zeros(length(t), N_execucoes); % ADICIONADO: Matriz para MSE de y_dot
mse_z     = zeros(length(t), N_execucoes);
mse_omega = zeros(length(t), N_execucoes);

% Vetores reais para comparação
y_real_vetor     = A_real * sin(omega_real * t);
ydot_real_vetor  = A_real * omega_real * cos(omega_real * t);
z_real_vetor     = z_real * ones(length(t), 1);
omega_real_vetor = omega_real * ones(length(t), 1);

for exec = 1:N_execucoes
    % Condições Iniciais do Filtro (com erro proposital)
    x_est = [10,10,-10]; % z começa em 1.0 [y; ydot; z]
    P = eye(grau) * 0.5;

    est_param = zeros(length(t), grau);
    omega_est_vetor = zeros(length(t), 1);

    for i = 1:length(t)
        % --- GERAÇÃO DA MEDIÇÃO REAL ---
        z_med = y_real_vetor(i) + sqrt(R) * randn();

        % Estados atuais para as Jacobianas
        y_hat  = x_est(1);
        yd_hat = x_est(2);
        z_hat  = x_est(3);

        % --- EKF: PREDIÇÃO (Modelo 5) ---
        x_pred = [ y_hat + yd_hat * dt;
                   yd_hat - (z_hat * y_hat) * dt;
                   z_hat ];

        % Fk e Qk adaptativos (Dependem de y_hat e z_hat)
        Fk = [ 1,        dt,    0;
              -z_hat*dt, 1,    -y_hat*dt;
               0,        0,     1 ];

        Qk = sig_z_sq * [ 0, 0, 0;
            0, (y_hat^2 * (dt^3)/3), (-y_hat * (dt^2)/2);
            0, (-y_hat * (dt^2)/2), dt ];

        P_pred = Fk * P * Fk' + Qk;

        % --- EKF: CORREÇÃO ---
        y_tilde = z_med - H * x_pred;
        S = H * P_pred * H' + R;
        K = (P_pred * H') / S;

        x_est = x_pred + K * y_tilde;
        P = (eye(grau) - K * H) * P_pred;

        % Salva resultados
        est_param(i, :) = x_est';
        % Usa abs para garantir estabilidade numérica na raiz
        omega_est_vetor(i) = sqrt(abs(x_est(3)));
    end

    % --- CÁLCULO DOS ERROS DA EXECUÇÃO ATUAL ---
    err_y     = y_real_vetor - est_param(:, 1);
    err_ydot  = ydot_real_vetor - est_param(:, 2); % ADICIONADO: Erro de y_dot
    err_z     = z_real_vetor - est_param(:, 3);
    err_omega = omega_real_vetor - omega_est_vetor;

    todos_mse(:, exec) = cumsum(err_y.^2) ./ (1:length(t))';
    mse_y(:, exec)     = todos_mse(:, exec);
    mse_ydot(:, exec)  = cumsum(err_ydot.^2) ./ (1:length(t))'; % ADICIONADO: MSE progressivo de y_dot
    mse_z(:, exec)     = cumsum(err_z.^2) ./ (1:length(t))';
    mse_omega(:, exec) = cumsum(err_omega.^2) ./ (1:length(t))';
end

% Médias das execuções
mse_medio_geral = mean(todos_mse, 2);
mse_medio_ydot  = mean(mse_ydot, 2); % ADICIONADO: Média do MSE de y_dot
mse_medio_z     = mean(mse_z, 2);
mse_medio_omega = mean(mse_omega, 2);

% Definição global de tamanho de fonte
fSize = 14;

% ---- 4. Plotagem ----

% Figura 1: Rastreamento e Convergência (Última Execução)
figure(1);
subplot(3,1,1);
plot(t, y_real_vetor, 'b', 'LineWidth', 2); hold on;
plot(t, est_param(:,1), 'r--', 'LineWidth', 1.5);
title('1. Rastreamento do Sinal ($y$)', 'FontSize', fSize, 'Interpreter', 'latex');
ylabel('Amplitude', 'FontSize', fSize);
legend({'Real','EKF'}, 'FontSize', fSize-2, 'Location', 'northeast');
grid on;
set(gca, 'FontSize', fSize);

subplot(3,1,2);
plot(t, ydot_real_vetor, 'k--', 'LineWidth', 2); hold on;
plot(t, est_param(:,2), 'g', 'LineWidth', 2);
title('2. Convergência de ($\dot{y}$)', 'FontSize', fSize, 'Interpreter', 'latex');
ylabel('Taxa', 'FontSize', fSize, 'Interpreter', 'latex');
grid on;
set(gca, 'FontSize', fSize);

subplot(3,1,3);
plot(t, omega_real_vetor, 'k--', 'LineWidth', 2); hold on;
plot(t, omega_est_vetor, 'm', 'LineWidth', 2);
title('3. Frequência Recuperada ($\omega$)', 'FontSize', fSize, 'Interpreter', 'latex');
ylabel('rad/s', 'FontSize', fSize);
xlabel('Tempo (s)', 'FontSize', fSize);
grid on;
set(gca, 'FontSize', fSize);

% Figura 2: Análise de Robustez - MSE do Sinal (Estilo Modelo 2)
figure(2);
plot(t, todos_mse, 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5); hold on;
plot(t, mse_medio_geral, 'k', 'LineWidth', 2.5);
title(sprintf('Erro Quadrático Médio (MSE) do Sinal - %d Execuções', N_execucoes), 'FontSize', fSize);
xlabel('Tempo (s)', 'FontSize', fSize);
ylabel('MSE Progressivo', 'FontSize', fSize);
legend({'Execuções Individuais', 'Média Geral'}, 'FontSize', fSize-2, 'Location', 'northeast');
grid on;
set(gca, 'FontSize', fSize);

% Figura 3: MSE Médio por Estado - ATUALIZADO PARA 3 SUBPLOTS
figure(3);
subplot(3,1,1); % MODIFICADO: De (2,1,1) para (3,1,1)
plot(t, mse_medio_ydot, 'm', 'LineWidth', 2.5); % ADICIONADO: Plot de y_dot
title('MSE Médio: Derivada do Sinal ($\dot{y}$)', 'FontSize', fSize, 'Interpreter', 'latex');
ylabel('MSE', 'FontSize', fSize);
grid on;
set(gca, 'FontSize', fSize);

subplot(3,1,2); % MODIFICADO: De (2,1,2) para (3,1,2)
plot(t, mse_medio_z, 'g', 'LineWidth', 2.5);
title('MSE Médio: Frequência ao Quadrado ($z$)', 'FontSize', fSize, 'Interpreter', 'latex');
ylabel('MSE', 'FontSize', fSize);
grid on;
set(gca, 'FontSize', fSize);

subplot(3,1,3); % ADICIONADO: Terceiro subplot
plot(t, mse_medio_omega, 'r', 'LineWidth', 2.5);
title('MSE Médio: Frequência Angular ($\omega$)', 'FontSize', fSize, 'Interpreter', 'latex');
ylabel('MSE', 'FontSize', fSize);
xlabel('Tempo (s)', 'FontSize', fSize);
grid on;
set(gca, 'FontSize', fSize);
