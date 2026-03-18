% Sine Wave Generation based on Six Models
clear; clc; close all;
pkg load control % Carrega o pacote de controle no Octave

% Common parameters
A_val = 2;          % Amplitude
f = 1;              % Frequency in Hz
omega_val = 2*pi*f; % Angular frequency
phi_val = pi/4;     % Phase (for Model 6)
tspan = [0 10];     % Time span for simulation (seconds)
t_eval = linspace(tspan(1), tspan(2), 500);

%% Model 1
% x = [phi; omega], phi_dot = omega, omega_dot = 0
% y = A * sin(phi)
f1 = @(t, x) [x(2); 0];
x0_1 = [0; omega_val];
[t1, x1] = ode45(f1, t_eval, x0_1);
y1 = A_val * sin(x1(:, 1));

%% Model 2
% x = [phi; omega; A], phi_dot = omega, omega_dot = 0, A_dot = 0
% y = A * sin(phi)
f2 = @(t, x) [x(2); 0; 0];
x0_2 = [0; omega_val; A_val];
[t2, x2] = ode45(f2, t_eval, x0_2);
y2 = x2(:, 3) .* sin(x2(:, 1));

%% Model 3
% x = [y; y_dot; omega], x_dot = [x2; -(x3^2)*x1; 0]
f3 = @(t, x) [x(2); -(x(3)^2)*x(1); 0];
x0_3 = [0; A_val * omega_val; omega_val]; % y(0)=0, y_dot(0)=A*omega
[t3, x3] = ode45(f3, t_eval, x0_3);
y3 = x3(:, 1);

%% Model 4
% x = [y; y_dot; mu], mu = omega^2, x_dot = [x2; -x3*x1; 0]
f4 = @(t, x) [x(2); -x(3)*x(1); 0];
x0_4 = [0; A_val * omega_val; omega_val^2];
[t4, x4] = ode45(f4, t_eval, x0_4);
y4 = x4(:, 1);

%% Model 5
% x_dot = Ax, A = [0 1; -omega^2 0], C = [1 0]
% Initial state x0 = [0; omega*A] for y = A*sin(omega*t)
A_mat = [0 1; -omega_val^2 0];
B_mat = [0; 0]; % GNU Octave exige a matriz declarada
C_mat = [1 0];  % Corrigido para [1 0] para não escalar a amplitude
D_mat = 0;      % GNU Octave exige a matriz declarada
x0_5 = [0; omega_val * A_val];
sys5 = ss(A_mat, B_mat, C_mat, D_mat);
[y5, t5] = initial(sys5, x0_5, t_eval);

%% Model 6
% Same A, B, C, D as Model 5, but different initial state
% x0 = [A*sin(phi); A*omega*cos(phi)]
x0_6 = [A_val * sin(phi_val); A_val * omega_val * cos(phi_val)];
sys6 = ss(A_mat, B_mat, C_mat, D_mat);
[y6, t6] = initial(sys6, x0_6, t_eval);

%% Plotting Results
figure('Name', 'Sine Wave Models Comparison');

subplot(3,2,1); plot(t1, y1, 'linewidth', 1.5); title('Model 1'); grid on;
subplot(3,2,2); plot(t2, y2, 'linewidth', 1.5); title('Model 2'); grid on;
subplot(3,2,3); plot(t3, y3, 'linewidth', 1.5); title('Model 3'); grid on;
subplot(3,2,4); plot(t4, y4, 'linewidth', 1.5); title('Model 4'); grid on;
subplot(3,2,5); plot(t5, y5, 'linewidth', 1.5); title('Model 5'); grid on;
subplot(3,2,6); plot(t6, y6, 'linewidth', 1.5); title('Model 6 (Phase Shift)'); grid on;
xlabel('Time (s)');
