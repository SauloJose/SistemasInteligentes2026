function [res_real, res_est, z_medido, mse_progressivo] = filtro_kalman_generico(f_func, h_func, F_jacob, H_jacob, Q, R, x0_real, x0_est, t)
    % filtro_kalman_generico: Implementação do Filtro de Kalman (Linear ou Estendido)
    % Entradas:
    %   f_func  - Função de transição de estado não-linear f(x)
    %   h_func  - Função de medição não-linear h(x)
    %   F_jacob - Função que retorna o Jacobiano de f avaliado em x
    %   H_jacob - Função que retorna o Jacobiano de h avaliado em x
    %   Q, R    - Matrizes de covariância de processo e medição
    %   x0_real - Estado inicial real do sistema
    %   x0_est  - Estimativa inicial do filtro
    %   t       - Vetor de tempo

    num_steps = length(t);
    n_states = length(x0_est);
    num_sensors = size(R, 1);

    % Inicialização
    P = eye(n_states);
    x_real = x0_real;
    x_est = x0_est;

    % Pré-alocação para velocidade computacional
    res_real = zeros(num_steps, n_states);
    res_est  = zeros(num_steps, n_states);
    z_medido = zeros(num_steps, num_sensors);

    for i = 1:num_steps
        % 1. Salva o estado real ATUAL
        res_real(i, :) = x_real';

        % 2. Gera a medição baseada no estado real (com ruído)
        ruido_sensor = sqrt(R) * randn(num_sensors, 1);
        z = h_func(x_real) + ruido_sensor;
        z_medido(i, :) = z';

        % --- ETAPA KALMAN: Predição ---
        x_est_prev = f_func(x_est);       % Propagação do estado
        F = F_jacob(x_est);               % Avalia o Jacobiano no ponto atual
        P_prev = F * P * F' + Q;          % Propagação da covariância

        % --- ETAPA KALMAN: Atualização/Correção ---
        H = H_jacob(x_est_prev);          % Avalia o Jacobiano da medição
        y_res = z - h_func(x_est_prev);   % Inovação (Erro residual)

        % Cálculo do Ganho K (Usando '/' no lugar de 'inv' para estabilidade numérica)
        S = H * P_prev * H' + R;
        K = P_prev * H' / S;

        x_est = x_est_prev + K * y_res;
        P = (eye(n_states) - K * H) * P_prev;

        % 5. Salva a estimativa corrigida
        res_est(i, :) = x_est';

        % 6. Evolui o sistema real (Planta)
        x_real = f_func(x_real);
    end

    % --- Cálculo do Erro Quadrático Médio (MSE) ---
    erro_puro = res_real(:,1) - res_est(:,1);
    erro_quadrado = erro_puro.^2;
    mse_progressivo = cumsum(erro_quadrado) ./ (1:num_steps)';

    % --- Plotagem Embutida ---
    figure('Name', 'Resultados do Filtro de Kalman Generico');

    subplot(3,1,1);
    plot(t, res_real(:,1), 'b', 'linewidth', 2); hold on;
    plot(t, z_medido(:,1), 'r.', 'markersize', 6); hold off;
    title('1. Sistema Real vs Medição Ruidosa');
    grid on;

    subplot(3,1,2);
    plot(t, res_real(:,1), 'b', 'linewidth', 2); hold on;
    plot(t, res_est(:,1), 'g-', 'linewidth', 1.5); hold off;
    title('2. Rastreamento do Filtro de Kalman');
    legend('Real', 'Estimado');
    grid on;

    subplot(3,1,3);
    plot(t, mse_progressivo, 'm', 'linewidth', 2);
    title('3. Erro Quadrático Médio (MSE) Progressivo');
    xlabel('Tempo (s)'); ylabel('MSE');
    grid on;
end
