import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import expm
from numba import njit
from typing import Tuple
from numpy.linalg import inv

# Configuração de estilo para os gráficos ficarem mais legíveis
plt.style.use('seaborn-v0_8-whitegrid')

#====================================================================
# Funções básicas
#====================================================================
def dist(point_a: np.ndarray, point_b: np.ndarray) -> float:
    """
    Calcula a distância euclidiana entre dois pontos n-dimensionais.
    """
    return np.linalg.norm(point_a - point_b)


def h_obs(pos: np.ndarray, bases: np.ndarray) -> np.ndarray:
    """
    Calcula o vetor de observação h(x) (distâncias teóricas).
    
    Retorna um vetor coluna (N, 1) contendo as distâncias euclidianas 
    entre a posição atual e as estações base.
    """
    # Cálculo vetorizado das distâncias (Norma L2 por linha)
    z = np.linalg.norm(bases - pos, axis=1)
    
    return z.reshape(-1, 1)


# Gerar valores estocásticos
def create_cov(sigmas: np.ndarray) -> np.ndarray:
    """
    Gera uma matriz de covariância diagonal a partir de uma lista de desvios padrões.
    Útil para inicializar matrizes de ruído independentes (ex: R ou Q diagonal).
    """
    # Transformamos os desvios em variâncias (quadrado) e criamos a diagonal 
    return np.diag(np.square(sigmas))

def get_noise_vec(R: np.ndarray) -> np.ndarray:
    """
    Gera uma amostra de vetor de ruído branco (média zero) baseado na matriz de covariância.
    Retorna um vetor coluna (N, 1).
    """
    # Determina o número de sensores (dimensão da matriz)
    n = R.shape[0]
    
    # Gera o ruído com média zero para todos os sensores
    v = np.random.multivariate_normal(np.zeros(n), R)
    
    # Retorna como vetor coluna (n, 1) para manter o padrão do filtro
    return v.reshape(-1, 1)

#====================================================================
# Funções do modelo PV
#====================================================================
@njit
def get_Q_disc_PV(q_diag: list, dt: float) -> np.ndarray:
    """
    q_diag: lista com as variâncias contínuas [q1(px), q2(py), q3(vx), q4(vy)]
    dt: delta t (tempo de amostragem)
    """
    q1, q2, q3, q4 = q_diag
    
    Q = np.zeros((4, 4))
    
    # Eixo X
    Q[0, 0] = q1 * dt + q3 * (dt**3) / 3.0
    Q[0, 2] = q3 * (dt**2) / 2.0
    Q[2, 0] = Q[0, 2]
    Q[2, 2] = q3 * dt
    
    # Eixo Y
    Q[1, 1] = q2 * dt + q4 * (dt**3) / 3.0
    Q[1, 3] = q4 * (dt**2) / 2.0
    Q[3, 1] = Q[1, 3]
    Q[3, 3] = q4 * dt
    
    return Q


@njit
def calc_h_and_H_PV(x: np.ndarray, bases: np.ndarray):
    """
    Calcula as distâncias esperadas h(x) e a Matriz Jacobiana H para o modelo PV.
    x: Vetor de estado (4, 1) -> [px, py, vx, vy]
    """
    n_sensores = bases.shape[0]
    hx = np.zeros((n_sensores, 1))
    H = np.zeros((n_sensores, 4)) # 4 colunas (px, py, vx, vy)
    
    px, py = x[0, 0], x[1, 0]
    
    for i in range(n_sensores):
        bx, by = bases[i, 0], bases[i, 1]
        dx = px - bx
        dy = py - by
        dist = np.sqrt(dx**2 + dy**2)
        
        hx[i, 0] = dist
        
        # Derivadas parciais (Jacobiana)
        # Evitamos divisão por zero caso o alvo esteja exatamente sobre a antena
        if dist > 1e-8:
            H[i, 0] = dx / dist  # d(dist)/d(px)
            H[i, 1] = dy / dist  # d(dist)/d(py)
        # H[i, 2] e H[i, 3] continuam 0 pois a distância não depende da velocidade

    return hx, H


class EKF_PV:
    def __init__(self, x0: np.ndarray, P0: np.ndarray, F: np.ndarray, Q: np.ndarray, R: np.ndarray):
        """
        Inicializa o Filtro de Kalman Estendido PV.
        """
        self.x = x0.copy()  # Vetor de estado (4x1)
        self.P = P0.copy()  # Covariância do erro (4x4)
        self.F = F          # Matriz de transição discreta (4x4)
        self.Q = Q          # Covariância do processo discreta (4x4)
        self.R = R          # Covariância da medição (NxN)
        self.I = np.eye(4)  # Matriz Identidade para o Update
        
    def predict(self):
        """ Passo 1: Predição do Estado e da Covariância """
        self.x = self.F @ self.x
        self.P = self.F @ self.P @ self.F.T + self.Q
        
    def update(self, z: np.ndarray, bases: np.ndarray):
        """ Passo 2: Correção baseada nas medições (z) """
        # Calcula medição esperada h(x) e Jacobiana (H) no estado previsto atual
        hx, H = calc_h_and_H_PV(self.x, bases)
        
        # Inovação (Residual)
        z = z.reshape(-1, 1) # Garante que z seja um vetor coluna
        y = z - hx
        
        # Covariância da Inovação (S) e Ganho de Kalman (K)
        S = H @ self.P @ H.T + self.R
        K = self.P @ H.T @ inv(S)
        
        # Atualização de Estado e Covariância
        self.x = self.x + K @ y
        self.P = (self.I - K @ H) @ self.P


#====================================================================
# Funções do modelo PVA
#====================================================================
@njit
def get_Q_disc_PVA(q_diag: np.ndarray, dt: float) -> np.ndarray:
    """
    Gera a matriz Q discreta (6x6) para o modelo PVA generalizado.
    q_diag: array com as variâncias contínuas [q1(px), q2(py), q3(vx), q4(vy), q5(ax), q6(ay)]
    dt: delta t (tempo de amostragem)
    """
    q1, q2, q3, q4, q5, q6 = q_diag
    
    Q = np.zeros((6, 6))
    
    # --- Eixo X (Índices 0, 2, 4 correspondentes a px, vx, ax) ---
    Q[0, 0] = q1 * dt + q3 * (dt**3) / 3.0 + q5 * (dt**5) / 20.0
    Q[0, 2] = q3 * (dt**2) / 2.0 + q5 * (dt**4) / 8.0
    Q[2, 0] = Q[0, 2]
    
    Q[0, 4] = q5 * (dt**3) / 6.0
    Q[4, 0] = Q[0, 4]
    
    Q[2, 2] = q3 * dt + q5 * (dt**3) / 3.0
    Q[2, 4] = q5 * (dt**2) / 2.0
    Q[4, 2] = Q[2, 4]
    
    Q[4, 4] = q5 * dt
    
    # --- Eixo Y (Índices 1, 3, 5 correspondentes a py, vy, ay) ---
    Q[1, 1] = q2 * dt + q4 * (dt**3) / 3.0 + q6 * (dt**5) / 20.0
    Q[1, 3] = q4 * (dt**2) / 2.0 + q6 * (dt**4) / 8.0
    Q[3, 1] = Q[1, 3]
    
    Q[1, 5] = q6 * (dt**3) / 6.0
    Q[5, 1] = Q[1, 5]
    
    Q[3, 3] = q4 * dt + q6 * (dt**3) / 3.0
    Q[3, 5] = q6 * (dt**2) / 2.0
    Q[5, 3] = Q[3, 5]
    
    Q[5, 5] = q6 * dt
    
    return Q

@njit
def calc_h_and_H_PVA(x: np.ndarray, bases: np.ndarray):
    """
    Calcula h(x) e a Matriz Jacobiana H para o modelo PVA.
    x: Vetor de estado (6, 1) -> [px, py, vx, vy, ax, ay]
    """
    n_sensores = bases.shape[0]
    hx = np.zeros((n_sensores, 1))
    H = np.zeros((n_sensores, 6)) # 6 colunas
    
    px, py = x[0, 0], x[1, 0]
    
    for i in range(n_sensores):
        bx, by = bases[i, 0], bases[i, 1]
        dx = px - bx
        dy = py - by
        dist = np.sqrt(dx**2 + dy**2)
        
        hx[i, 0] = dist
        
        if dist > 1e-8:
            H[i, 0] = dx / dist
            H[i, 1] = dy / dist
        # Velocidades e acelerações não afetam diretamente a medição da distância, ficam em 0.

    return hx, H

class EKF_PVA:
    def __init__(self, x0: np.ndarray, P0: np.ndarray, F: np.ndarray, Q: np.ndarray, R: np.ndarray):
        """
        Inicializa o Filtro de Kalman Estendido PVA.
        """
        self.x = x0.copy()  # Vetor de estado (6x1)
        self.P = P0.copy()  # Covariância do erro (6x6)
        self.F = F          # Matriz de transição discreta (6x6)
        self.Q = Q          # Covariância do processo (6x6)
        self.R = R          # Covariância da medição (NxN)
        self.I = np.eye(6)  # Matriz Identidade para o Update
        
    def predict(self):
        """ Passo 1: Predição do Estado e da Covariância """
        self.x = self.F @ self.x
        self.P = self.F @ self.P @ self.F.T + self.Q
        
    def update(self, z: np.ndarray, bases: np.ndarray):
        """ Passo 2: Correção baseada nas medições (z) """
        hx, H = calc_h_and_H_PVA(self.x, bases)

        z = z.reshape(-1, 1) # Garante que z seja um vetor coluna
        y = z - hx
        
        S = H @ self.P @ H.T + self.R
        K = self.P @ H.T @ inv(S)
        
        self.x = self.x + K @ y
        self.P = (self.I - K @ H) @ self.P

#====================================================================
# Gerador de trajetórias para utilizar o Filtro de Kalman
#====================================================================
class TrajectoryGenerator:
    def __init__(self, dt: float, bases: np.ndarray):
        """
        Inicializa o gerador de trajetórias.
        dt: Intervalo de tempo entre amostras (s)
        bases: Posição das antenas base (N, 2) [[bx1, by1], ...]
        """
        self.dt = dt
        self.bases = bases
        
        # Plot styling
        plt.style.use('seaborn-v0_8-whitegrid')
        self.fig_size = (10, 8)

    def _init_vectors(self, num_steps: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Inicializa os vetores de tempo e estado PVA (6D)."""
        t = np.linspace(0, (num_steps - 1) * self.dt, num_steps)
        # Estados Verdadeiros (Ground Truth): [px, py, vx, vy, ax, ay]
        states_gt = np.zeros((num_steps, 6))
        return t, states_gt

    def generate_circle(self, radius: float, center: Tuple[float, float], 
                        linear_velocity: float, duration: float) -> Tuple[np.ndarray, np.ndarray]:
        """
        Gera uma trajetória circular com cinemática PVA completa.
        A aceleração gerada é a aceleração centrípeta.
        """
        num_steps = int(duration / self.dt)
        t, states = self._init_vectors(num_steps)
        cx, cy = center
        
        # Velocidade angular (w = v / r)
        omega = linear_velocity / radius
        
        # Ângulo em função do tempo (theta = w * t)
        theta = omega * t
        
        # Posições: [px, py]
        states[:, 0] = cx + radius * np.cos(theta)
        states[:, 1] = cy + radius * np.sin(theta)
        
        # Velocidades (derivada da posição): [vx, vy]
        states[:, 2] = -radius * omega * np.sin(theta)
        states[:, 3] =  radius * omega * np.cos(theta)
        
        # Acelerações (derivada da velocidade): [ax, ay]
        states[:, 4] = -radius * (omega**2) * np.cos(theta)
        states[:, 5] = -radius * (omega**2) * np.sin(theta)
        
        return t, states

    def generate_square(self, side_length: float, bottom_left: Tuple[float, float], 
                        linear_velocity: float) -> Tuple[np.ndarray, np.ndarray]:
        """
        Gera uma trajetória quadrada. 
        Nota: Nas quinas, há descontinuidade de velocidade (aceleração infinita teórica).
        Aqui simulamos velocidade constante nas arestas e mudança brusca na quina.
        """
        time_per_side = side_length / linear_velocity
        steps_per_side = int(time_per_side / self.dt)
        num_steps = steps_per_side * 4
        
        t, states = self._init_vectors(num_steps)
        x0, y0 = bottom_left
        v = linear_velocity
        
        # Areste 1: Direita (v_x = v, v_y = 0)
        idx = slice(0, steps_per_side)
        states[idx, 0] = x0 + v * t[idx]
        states[idx, 1] = y0
        states[idx, 2] = v
        
        # Areste 2: Cima (v_x = 0, v_y = v)
        idx = slice(steps_per_side, 2 * steps_per_side)
        states[idx, 0] = x0 + side_length
        states[idx, 1] = y0 + v * (t[idx] - time_per_side)
        states[idx, 3] = v
        
        # Areste 3: Esquerda (v_x = -v, v_y = 0)
        idx = slice(2 * steps_per_side, 3 * steps_per_side)
        states[idx, 0] = x0 + side_length - v * (t[idx] - 2 * time_per_side)
        states[idx, 1] = y0 + side_length
        states[idx, 2] = -v
        
        # Areste 4: Baixo (v_x = 0, v_y = -v)
        idx = slice(3 * steps_per_side, 4 * steps_per_side)
        states[idx, 0] = x0
        states[idx, 1] = y0 + side_length - v * (t[idx] - 3 * time_per_side)
        states[idx, 3] = -v
        
        # Acelerações são zero nas arestas (modelo PV perfeito).
        # As quinas geram picos de erro no filtro PV.
        return t, states

    def generate_tanh_curve(self, start_pos: Tuple[float, float], end_pos: Tuple[float, float], 
                            amplitude_y: float, smoothness: float, duration: float) -> Tuple[np.ndarray, np.ndarray]:
        """
        Gera uma trajetória suave baseada na tangente hiperbólica (curva em S),
        semelhante a manobras de troca de faixa, com PVA completo.
        """
        num_steps = int(duration / self.dt)
        t, states = self._init_vectors(num_steps)
        x_start, y_start = start_pos
        x_end, y_end = end_pos
        
        # Tempo normalizado (0 a 1) e centralizado (-0.5 a 0.5)
        t_norm = (t / duration) - 0.5
        
        # --- Eixo X (Movimento Linear constante) ---
        dist_x = x_end - x_start
        v_x_const = dist_x / duration
        states[:, 0] = x_start + v_x_const * t # Pos X
        states[:, 2] = v_x_const                # Vel X
        # Acc X = 0
        # --- Eixo Y (Curva Tanh) ---
        y_mid = (y_end + y_start) / 2.0
        
        # Argumento da tanh: controla a inclinação da curva
        # arg = smoothness * t_norm
        
        # Posição Y: y = y_mid + Amp * tanh(smoothness * t_norm)
        tanh_arg = smoothness * t_norm
        states[:, 1] = y_mid + amplitude_y * np.tanh(tanh_arg)
        
        # Derivadas analíticas para Velocidade e Aceleração em Y
        # d/dx tanh(x) = sech^2(x) = 1 - tanh^2(x)
        # Necessário usar regra da cadeia: d(smoothness*t_norm)/dt = smoothness / duration
        sech2_val = 1.0 - (np.tanh(tanh_arg))**2
        factor_dt = smoothness / duration
        
        # Velocidade Y: v_y = Amp * factor_dt * sech^2(arg)
        states[:, 3] = amplitude_y * factor_dt * sech2_val
        
        # Aceleração Y: a_y = d/dt v_y = Amp * factor_dt^2 * d/d(arg)sech^2(arg)
        # d/dx sech^2(x) = -2 * sech^2(x) * tanh(x)
        states[:, 5] = amplitude_y * (factor_dt**2) * (-2.0 * sech2_val * np.tanh(tanh_arg))
        
        return t, states

    def plot_scenario(self, t: np.ndarray, states_gt: np.ndarray, title: str):
        """Planta a trajetória (Ground Truth) e as bases."""
        px = states_gt[:, 0]
        py = states_gt[:, 1]
        vx = states_gt[:, 2]
        vy = states_gt[:, 3]
        ax = states_gt[:, 4]
        ay = states_gt[:, 5]
        
        # Cria figura com subplots (Trajetória + Cinematica)
        fig = plt.figure(figsize=(14, 10))
        gs = fig.add_gridspec(3, 2)
        
        # --- Plot 1: Espacial (X vs Y) ---
        ax_space = fig.add_subplot(gs[0:2, 0])
        
        # Plot Trajetória ideal
        ax_space.plot(px, py, 'g-', linewidth=2.5, label='Trajetória Ideal (GT)')
        
        # Início e Fim
        ax_space.plot(px[0], py[0], 'go', markersize=8, label='Início')
        ax_space.plot(px[-1], py[-1], 'gx', markersize=10, mew=2, label='Fim')
        
        # Plot Bases
        ax_space.scatter(self.bases[:, 0], self.bases[:, 1], 
                         marker='^', s=150, c='red', edgecolor='black', 
                         label='Antenas Base', zorder=5)
        # Numeração das bases
        for i, (bx, by) in enumerate(self.bases):
            ax_space.annotate(f'B{i}', (bx, by), xytext=(5, 5), textcoords='offset points', fontsize=12, fontweight='bold')
            
        ax_space.set_title(f'Cenário Espacial: {title}', fontsize=14)
        ax_space.set_xlabel('Posição X (m)', fontsize=12)
        ax_space.set_ylabel('Posição Y (m)', fontsize=12)
        ax_space.legend(loc='best', frameon=True, shadow=True)
        ax_space.axis('equal') # Mantém proporção 1:1 m
        ax_space.grid(True, which='both', linestyle='--')

        # --- Subplots Cinemática (Tempo) ---
        # Velocidade
        ax_vel = fig.add_subplot(gs[0, 1])
        ax_vel.plot(t, vx, 'b--', label='$v_x$')
        ax_vel.plot(t, vy, 'r--', label='$v_y$')
        ax_vel.set_title('Velocidade Ideal', fontsize=12)
        ax_vel.set_ylabel('(m/s)')
        ax_vel.legend()
        ax_vel.grid(True)
        
        # Aceleração
        ax_acc = fig.add_subplot(gs[1, 1])
        ax_acc.plot(t, ax, 'b-', label='$a_x$')
        ax_acc.plot(t, ay, 'r-', label='$a_y$')
        ax_acc.set_title('Aceleração Ideal (Necessária para a curva)', fontsize=12)
        ax_acc.set_ylabel('(m/$s^2$)')
        ax_acc.legend()
        ax_acc.grid(True)
        
        # Tanh Y específica (para mostrar a forma)
        ax_tanh = fig.add_subplot(gs[2, 1])
        ax_tanh.plot(t, py, 'k-', label='$p_y$ (Forma Tanh)')
        ax_tanh.set_title('Perfil de Posição Y (Forma da Tanh)', fontsize=12)
        ax_tanh.set_xlabel('Tempo (s)')
        ax_tanh.set_ylabel('(m)')
        ax_tanh.grid(True)
        
        fig.tight_layout()
        plt.show()

# ====================================================================
# EXEMPLO DE USO DO MODELO, nesse caso PVA
# ====================================================================

if __name__ == "__main__":
    #===========================================
    # Definindo parâmetros globais da simulação
    #===========================================
    dT_sim = 1e-1
    
    # Posição das estações base (X,Y).
    S1 = np.array([0,0]).T
    S2 = np.array([600,600]).T
    S3 = np.array([0,600]).T
    S4 = np.array([600,0]).T
    
    #Variável para iterar sobre elas
    bases = np.array([S1,S2,S3,S4])
    
    # Número de antenas para medição
    N = bases.shape[0]
    
    print("Tem",N, "Antenas")
    
    antenas_base = bases
    gen = TrajectoryGenerator(dT_sim, antenas_base) 

    #===========================================
    # Gerando Trajetórias
    #===========================================
    t_sq, pva_sq = gen.generate_square(side_length=400, bottom_left=(100, 100), linear_velocity=8)
    
    # Puxando valores reais 
    # pva_sq possui 6 colunas, mas para o PV usamos apenas as 4 primeiras (px, py, vx, vy)
    px_gt = pva_sq[:, 0]
    py_gt = pva_sq[:, 1]
    vx_gt = pva_sq[:, 2]
    vy_gt = pva_sq[:, 3]
    
    #===========================================
    # Configurações do filtro de Kalman
    #===========================================
    # Estado inicial (assumimos que começa na posição correta, mas velocidade zerada)
    x0 = np.array([[px_gt[0]], [py_gt[0]], [vx_gt[0]], [vy_gt[0]]])
    P0 = np.eye(4) * 10.0  # Incerteza inicial
    
    # Matriz de Transição de Estado (F) para modelo PV
    F = np.array([
        [1, 0, dT_sim, 0],
        [0, 1, 0, dT_sim],
        [0, 0, 1,      0],
        [0, 0, 0,      1]
    ])
    
    # Covariância do Ruído de Processo (Q)
    # Variâncias do ruído contínuo
    q1 = 1
    q2 = 1
    q3 = 1
    q4 = 1
    
    # Criando matriz Q discreta:
    Q_disc = get_Q_disc_PV([q1, q2, q3, q4], dT_sim)
    
    # Covariância do Ruído de Medição (R)
    dp = 2  # seu valor ou objeto
    vetor = np.full(N, dp)
    
    # Pegando a matriz R (passando apenas os desvios)
    R = create_cov(vetor)
    
    # Instanciando o seu filtro
    ekf_pv = EKF_PV(x0, P0, F, Q_disc, R)
    
    # ==============================================
    # Histórico para salvar as estimativas
    num_steps = len(t_sq)
    estados_estimados = np.zeros((num_steps, 4))

    # ====================================================================
    # 3. LOOP DE SIMULAÇÃO E RASTREAMENTO
    # ====================================================================
    for k in range(num_steps):
        # a) Simulação da medição real (Extraímos a distância ideal e adicionamos ruído)
        # Criamos um vetor coluna com o estado real no instante k
        x_real_k = np.array([[px_gt[k]], [py_gt[k]], [vx_gt[k]], [vy_gt[k]]])
        pos = np.array([[px_gt[k]], [py_gt[k]]])
        
        # Usamos a sua função matemática para calcular as distâncias exatas h(x)
        hx_real, _ = calc_h_and_H_PV(x_real_k, antenas_base)
        
        # Injetamos ruído para simular o sensor Z real
        # Recuperando o vetor de erros aleatórios (4x1)
        v_t = get_noise_vec(R)
        
        # Medição final ruidosa (transpondo pos para broadcasting)
        z_medido = h_obs(pos.T, bases) + v_t
        
        # b) Execução do Filtro de Kalman
        ekf_pv.predict()
        ekf_pv.update(z_medido, antenas_base)
        
        # c) Salvando o estado estimado (achatando de (4,1) para (4,))
        estados_estimados[k, :] = ekf_pv.x.flatten()

    # Extraindo as estimativas para variáveis isoladas
    px_est = estados_estimados[:, 0]
    py_est = estados_estimados[:, 1]
    vx_est = estados_estimados[:, 2]
    vy_est = estados_estimados[:, 3]
    
    # ====================================================================
    # 4. PLOTAGEM E CÁLCULO DE ERROS
    # ====================================================================
    # Cálculo dos erros (Ground Truth - Estimado)
    erro_px = px_gt - px_est
    erro_py = py_gt - py_est
    erro_vx = vx_gt - vx_est
    erro_vy = vy_gt - vy_est
    
    # Plot 1: Trajetória 2D
    plt.figure(figsize=(10, 8))
    plt.plot(px_gt, py_gt, label='Trajetória Ideal (GT)', color='green', linewidth=2)
    plt.plot(px_est, py_est, label='EKF Estimativa (PV)', color='red', linestyle='dashed')
    plt.scatter(antenas_base[:, 0], antenas_base[:, 1], c='blue', marker='^', s=150, label='Antenas Base')
    for i, (bx, by) in enumerate(antenas_base):
        plt.annotate(f'B{i}', (bx, by), xytext=(5, 5), textcoords='offset points', fontsize=12, fontweight='bold')
    
    plt.title('Rastreamento EKF - Trajetória Quadrada (Modelo PV)', fontsize=14)
    plt.xlabel('Posição X (m)', fontsize=12)
    plt.ylabel('Posição Y (m)', fontsize=12)
    plt.legend()
    plt.grid(True, linestyle='--')
    plt.axis('equal')
    plt.show()
    
    # Plot 2: Análise de Erros nos Estados
    fig, axs = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Análise de Erro do EKF (Valor Ideal - Estimado)', fontsize=16)
    
    # Posição X
    axs[0, 0].plot(t_sq, erro_px, color='red')
    axs[0, 0].set_title('Erro em Posição X', fontsize=12)
    axs[0, 0].set_ylabel('Metros (m)')
    axs[0, 0].grid(True)
    
    # Posição Y
    axs[0, 1].plot(t_sq, erro_py, color='red')
    axs[0, 1].set_title('Erro em Posição Y', fontsize=12)
    axs[0, 1].grid(True)
    
    # Velocidade X
    axs[1, 0].plot(t_sq, erro_vx, color='blue')
    axs[1, 0].set_title('Erro em Velocidade X', fontsize=12)
    axs[1, 0].set_xlabel('Tempo (s)')
    axs[1, 0].set_ylabel('Velocidade (m/s)')
    axs[1, 0].grid(True)
    
    # Velocidade Y
    axs[1, 1].plot(t_sq, erro_vy, color='blue')
    axs[1, 1].set_title('Erro em Velocidade Y', fontsize=12)
    axs[1, 1].set_xlabel('Tempo (s)')
    axs[1, 1].grid(True)
    
    plt.tight_layout()
    plt.show()
    
    
    # RMS
    erro = erro_px + erro_py
    erro2 = np.square(erro)
    rms = np.mean(erro2)
    
    print("Erro médio quadrático: ", rms)

