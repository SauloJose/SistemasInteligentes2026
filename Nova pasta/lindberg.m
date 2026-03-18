clc
clear all
pkg load control % Carrega o pacote de controle no Octave


% 1. Sistema real (Velocidade de um carro)
% X é um oscilador harmônico.
omega = 2;
A = [0 1; -omega^2 0]; % Matriz que gera a senoide
B = [0; 0];            % Não precisa de entrada u, ele oscila sozinho
C = [1 0];             % Medimos a amplitude (x1)


% 2. Projeto do observador.
k1 = -2+3j;
k2 = k1';
polos_desejados = [k1;k2] %Criticamente amortecida
L = place(A', C', polos_desejados)';

auto_valores = eig(A-L*C)

% Simulação
dt = 0.01;
t = 0:dt:10;
u = ones(size(t));

% Estados iniciais e estimado
x_real = [1;0];
x_est = [0; 0];

% Valores puxados
res_real = [];
res_est = [];


% Looping no tempo para computar os dados
for i = 1:length(t)
  % Guardar dados para o gráfico
  res_real = [res_real; x_real'];
  res_est  = [res_est; x_est'];


  % Mundo real
  y = C*x_real;
  dx_real = A*x_real + B*u(i);
  x_real = x_real + dx_real*dt; %Integração numérica;

  % Observador:
  erro_leitura = y - (C*x_est);
  dx_est = A*x_est +B*u(i)+L*erro_leitura;
  x_est = x_est + dx_est*dt;

end

% 4. Erro total
erro_total = res_real - res_est;

% 5. PLOTAGEM
% Subplot de Estados
subplot(4,1,1);
plot(t, res_real(:,1), 'b', t, res_est(:,1), 'r--');
title('Posição (x1)'); legend('Real','Estimado'); grid on;

subplot(4,1,2);
plot(t, res_real(:,2), 'b', t, res_est(:,2), 'r--');
title('Velocidade (x2)'); legend('Real','Estimado'); grid on;

% Nova Figura para os Erros
subplot(4,1,3);
plot(t, erro_total(:,1), 'k', 'LineWidth', 1.5);
title('Erro de Estimação: Posição (e1)'); grid on;

subplot(4,1,4);
plot(t, erro_total(:,2), 'k', 'LineWidth', 1.5);
title('Erro de Estimação: Velocidade (e2)'); grid on;

fprintf('Convergência concluída. Erro final de posição: %.4f\n', erro_total(end,1));

