#============================================================
#
#
#
#
#
#============================================================

#Geração da matriz de Identicabilidade:
def minimal_polynomial_coeffs(A, tol=1e-8):
    """
    Coeficientes do polinômio mínimo de A: [a0, a1, ..., am], com a0 = 1.
 
    Resolve: A^m + a1·A^{m-1} + ... + am·I = 0
    Rearranjando: [vec(A^{m-1}) | ... | vec(I)] · [a1; ...; am] = -vec(A^m)
    """
    n = A.shape[0]
    powers = [np.eye(n)]
    for _ in range(n):
        powers.append(powers[-1] @ A)
 
    vec_powers = [p.flatten('F') for p in powers]
 
    for m in range(1, n + 1):
        # Colunas na ordem decrescente de potência: A^{m-1}, ..., I
        cols  = [vec_powers[j] for j in range(m - 1, -1, -1)]
        X_ls  = np.column_stack(cols)
        y_ls  = -vec_powers[m]
        coeffs, _, _, _ = np.linalg.lstsq(X_ls, y_ls, rcond=None)
 
        # BUG-2: usar norma explícita em vez de residuals[0]
        if np.linalg.norm(X_ls @ coeffs - y_ls) < tol:
            full = np.zeros(m + 1)
            full[0]  = 1.0
            full[1:] = coeffs   # BUG-1: sem [::-1]
            return np.real(full)
 
    return np.real(np.poly(A))   # fallback: polinômio característico

def compute_B_G_matrices(F, H, Gamma, W):
    """
    Calcula as listas B_list (i=1..m) e G_list (i=0..m) conforme definido no artigo.
    
    B_l = H · S_l · Γ,  S_l = Σ_{i=0}^{l-1} a_i F̄^{l-i-1}
    G_l = a_l·I − H · S_l · F · W,  G_0 = I
    
    Parâmetros:
    -----------
    F : ndarray (n_x, n_x)
    H : ndarray (n_z, n_x)
    Gamma : ndarray (n_x, n_v)
    W : ndarray (n_x, n_z)  # ganho do filtro (pode ser subótimo ou inicial)
    
    Retorna:
    --------
    B_list : list de ndarrays (comprimento m), cada um (n_z, n_v)
    G_list : list de ndarrays (comprimento m+1), cada um (n_z, n_z)
    a_coeffs : ndarray (m+1,) coeficientes do polinômio mínimo
    m : int
    """
    n_x  = F.shape[0]
    I_nx = np.eye(n_x)
    I_nz = np.eye(H.shape[0])
    F_bar = F @ (I_nx - W @ H)
 
    a = minimal_polynomial_coeffs(F_bar)
    a = a / a[0]        # garante a0 = 1
    m = len(a) - 1
 
    Fb_pow = [I_nx]
    for _ in range(1, m):
        Fb_pow.append(Fb_pow[-1] @ F_bar)
 
    S = [None] * (m + 1)
    for l in range(1, m + 1):
        S[l] = sum(a[i] * Fb_pow[l - i - 1] for i in range(l))
 
    B_list = [H @ S[l] @ Gamma                             for l in range(1, m + 1)]
    G_list = [I_nz] + [a[l] * I_nz - H @ S[l] @ F @ W    for l in range(1, m + 1)]
 
    return B_list, G_list, a, m
    
def build_identifiability_matrix(B_list, G_list, m, n_v, n_z):
    """
    Constrói a matriz de identificabilidade de covariância de ruído F
    conforme Algorithm 1 do artigo.

    Parâmetros:
    -----------
    B_list : list of ndarray
        Lista de matrizes B_i (i = 1..m), cada uma de dimensão (n_z, n_v).
    G_list : list of ndarray
        Lista de matrizes G_i (i = 0..m), cada uma de dimensão (n_z, n_z).
    m : int
        Ordem do polinômio mínimo da matriz F_bar = F*(I - W*H).
    n_v : int
        Dimensão do ruído de processo v(k).
    n_z : int
        Dimensão da medição z(k).

    Retorna:
    --------
    F_mat : ndarray
        Matriz de identificabilidade com dimensão ((m+1)*n_z^2, n_cols),
        onde n_cols = n_v*(n_v+1)//2 + n_z*(n_z+1)//2.
    """
    n_cols = n_v * (n_v + 1) // 2 + n_z * (n_z + 1) // 2
    I_mat  = np.zeros(((m + 1) * n_z ** 2, n_cols))
 
    Bc = [[B_list[i-1][:, c] for c in range(n_v)] for i in range(1, m+1)]
    Gc = [[G_list[i][:, c]   for c in range(n_z)] for i in range(m+1)]
 
    for j in range(m + 1):
        r, k = j * n_z**2, 0
        z2 = n_z**2
 
        for l in range(1, n_v + 1):
            b = (sum(np.outer(Bc[i-1][l-1], Bc[i-j-1][l-1])
                     for i in range(j+1, m+1))
                 if j < m else np.zeros((n_z, n_z)))
            I_mat[r:r+z2, k] = b.flatten('F');  k += 1
            for p in range(l+1, n_v+1):
                d = (sum(np.outer(Bc[i-1][l-1], Bc[i-j-1][p-1])
                         + np.outer(Bc[i-1][p-1], Bc[i-j-1][l-1])
                         for i in range(j+1, m+1))
                     if j < m else np.zeros((n_z, n_z)))
                I_mat[r:r+z2, k] = d.flatten('F');  k += 1
 
        for l in range(1, n_z + 1):
            g = sum(np.outer(Gc[i][l-1], Gc[i-j][l-1]) for i in range(j, m+1))
            I_mat[r:r+z2, k] = g.flatten('F');  k += 1
            for p in range(l+1, n_z+1):
                e = sum(np.outer(Gc[i][l-1], Gc[i-j][p-1])
                        + np.outer(Gc[i][p-1], Gc[i-j][l-1])
                        for i in range(j, m+1))
                I_mat[r:r+z2, k] = e.flatten('F');  k += 1
 
    return I_mat

#UTILS
def check_identifiability(F, H, Gamma, W=None, diagonal_Q=False, verbose=True):
    """Verifica identificabilidade de Q e R. W=None → usa W=0."""
    nx, nz, nv = F.shape[0], H.shape[0], Gamma.shape[1]
    if W is None:
        W = np.zeros((nx, nz))
 
    B, G, a, m = compute_B_G_matrices(F, H, Gamma, W)
    I_mat = build_identifiability_matrix(B, G, m, nv, nz)
    rank  = np.linalg.matrix_rank(I_mat)
    n_Q   = nv if diagonal_Q else nv * (nv + 1) // 2
    n_R   = nz * (nz + 1) // 2
 
    if verbose:
        print(f"Grau do polinômio mínimo m = {m}")
        print(f"Coeficientes a_i = {np.round(a, 4)}")
        print(f"Dimensão de I: {I_mat.shape}, rank = {rank}")
        print(f"Incógnitas: n_Q={n_Q}, n_R={n_R}, total={n_Q+n_R}")
        print(f"Identificável: {'SIM' if rank >= n_Q+n_R else 'NÃO'}")
 
    return rank >= n_Q + n_R, rank, I_mat
 
 
def normalized_innovation_squared(nu, S):
    """NIS médio: E[ν(k)^T S^{-1} ν(k)]  (Eq. 191)."""
    S_inv = inv(S)
    return float(np.mean([n @ S_inv @ n for n in nu]))

def init_gain_W(F,H,Gamma, Q0=None, R0=None):
    """
    Calcula um ganho de Kalman estabilizante inicial.
    
    X[k+1] = F @ X[k] + Gamma @ V[k] 
    Z[k]   = H @ X[k] + w[k] 
    
    Parâmetros:
    F, H, Gamma : Matrizes do Sistema (Já discretizadas)
    Q0, R0      : Estimativas iniciais das covariâncias 

    Retorna:
    W0: Ganho de Kalman estabilizante 
    """
    nx, nz, nv = F.shape[0], H.shape[0], Gamma.shape[1]
    Q0 = np.eye(nv) if Q0 is None else Q0
    R0 = np.eye(nz) if R0 is None else R0
 
    Q_eff = Gamma @ Q0 @ Gamma.T
    P_bar = solve_discrete_are(F.T, H.T, Q_eff, R0)
    S     = H @ P_bar @ H.T + R0
    W0    = P_bar @ H.T @ inv(S)
 
    if np.max(np.abs(eigvals(F @ (np.eye(nx) - W0 @ H)))) >= 1.0:
        print("Aviso: ganho inicial não estabilizante.")
    return W0


def generate_residuals(Z: np.ndarray, 
                       F: np.ndarray, 
                       H: np.ndarray, 
                       W: np.ndarray, 
                       x0_hat: np.ndarray = None):
    """
    Roda um Filtro de Kalman estacionário sobre os dados para gerar resíduos.
    
    Z : Array de medições com shape (N, nz) onde N é o número de amostras
    F, H : Matrizes do sistema
    W : Ganho de Kalman estabilizante
    x0_hat : Estimativa do estado inicial. Se None, assume zero.
    
    Retorna:
    v : Sequência de inovação (N, nz)
    mu : Sequência de resíduos pós-ajuste (N, nz)
    """
    N, nz  = Z.shape
    nx     = F.shape[0]
    nu     = np.zeros((N, nz))
    mu     = np.zeros((N, nz))
    x_hat  = np.zeros((nx, 1)) if x0_hat is None else x0_hat.reshape(nx, 1)
 
    for k in range(N):
        z_k   = Z[k].reshape(nz, 1)
        nu_k  = z_k - H @ x_hat
        x_upd = x_hat + W @ nu_k
        mu_k  = z_k - H @ x_upd
        nu[k] = nu_k.flatten()
        mu[k] = mu_k.flatten()
        x_hat = F @ x_upd
    return nu, mu



def compute_autocovariance(nu, M):
    """
    Ĉ(i) = 1/(N-M) · Σ_{j=1}^{N-M} ν(j+i) ν(j)^T,  i = 0 … M-1.
 
    BUG-4: divisor unificado N-M para todos os lags (conforme Eq. 50).
    """
    N, nz  = nu.shape
    denom  = max(N - M, 1)
    C      = np.zeros((M, nz, nz))
    for i in range(M):
        C[i] = nu[i: N - M + i].T @ nu[: N - M] / denom
    return C



#Definindo função objetivo
def objective_J(W, F, H, C):
    """
    Calcula a função objetivo J(W) conforme a Equação (54):
    J = ½ tr( Σ_{i=1}^{M-1} Θ(i) · X · E² · X^T )
    Θ(i) = Φ(i)^T E² Φ(i),  Φ(i) = H F̄^{i-1} F,  E² = diag(Ĉ(0))^{-1}

    Parâmetros
    ----------
    W : ndarray (nx, nz)
        Ganho do filtro de Kalman.
    F : ndarray (nx, nx)
        Matriz de transição de estado.
    H : ndarray (nz, nx)
        Matriz de observação.
    C : list de ndarray (M, nz, nz)
        Autocovariâncias amostrais das inovações.

    Retorna
    -------
    J : float
        Valor da função objetivo.
    """
    M       = len(C)
    nz, nx  = H.shape
    I_nx    = np.eye(nx)
    C0      = C[0]
    E2      = np.diag(1.0 / np.maximum(np.diag(C0), 1e-12))
    F_bar   = F @ (I_nx - W @ H)
    X       = compute_X(F_bar, F, H, C, M)
 
    J = 0.0
    for i in range(1, M):
        Phi_i   = H @ np.linalg.matrix_power(F_bar, i-1) @ F
        Theta_i = Phi_i.T @ E2 @ Phi_i
        J      += np.trace(Theta_i @ X @ E2 @ X.T)
    return 0.5 * J


def lstsq_reg(A, B, reg=1e-6):
    """Mínimos quadrados com regularização de Tikhonov."""
    return np.linalg.solve(A.T @ A + reg * np.eye(A.shape[1]), A.T @ B)

# Cálculo do X via pseudo-inversa
def compute_X(F_bar, F, H, C, M):
    A = np.vstack([H @ np.linalg.matrix_power(F_bar, i-1) @ F for i in range(1, M)])
    B = np.vstack([C[i] for i in range(1, M)])
    X, _, _, _ = np.linalg.lstsq(A, B, rcond=None)
    return X

# Calculo de Z via Lyapunov
def compute_Z(F_bar, Phi, E2, H, C, start, M):
    """
    Calcula Z para um dado índice inicial 'start', conforme a equação:
    Z = F_bar^T Z F_bar + Q_Z
    onde Q_Z = 1/2 * Σ_{j=start}^{M-1} [ (Phi(j)^T E^2 C(j) E^2 H) + (Phi(j)^T E^2 C(j) E^2 H)^T ]

    Parâmetros
    ----------
    F_bar : ndarray (nx, nx)
        Matriz de malha fechada.
    Phi : list de ndarray (M-1 elementos, cada (nz, nx))
        Lista pré-calculada de Phi(i) = H * F_bar^{i-1} * F.
    E2 : ndarray (nz, nz)
        Matriz de normalização E^2.
    H : ndarray (nz, nx)
        Matriz de observação.
    C : list de ndarray (M, nz, nz)
        Autocovariâncias amostrais.
    start : int
        Índice inicial do somatório (1 ≤ start ≤ M-1).
    M : int
        Número total de lags (len(C)).

    Retorna
    -------
    Z : ndarray (nx, nx)
        Solução da equação de Lyapunov.
    """
    nx = F_bar.shape[0]
    Q_Z = np.zeros((nx, nx))
    for j in range(start, M):          # j = start, ..., M-1
        term = Phi[j-1].T @ E2 @ C[j] @ E2 @ H   # (nx, nx)
        Q_Z += 0.5 * (term + term.T)
    Z = solve_discrete_lyapunov(F_bar.T, Q_Z)
    return Z

# Implementando função gradiente da função objetivo
def gradient_J(W, F, H, C):
    """
    Calcula o gradiente de J em relação a W.

        Parâmetros
    ----------
    W : numpy.ndarray, shape (nx, nz)
        Ganho estacionário do filtro de Kalman (a variável de otimização).
        nx = dimensão do estado, nz = dimensão da medição.

    F : numpy.ndarray, shape (nx, nx)
        Matriz de transição de estado do sistema (conhecida e fixa).

    H : numpy.ndarray, shape (nz, nx)
        Matriz de observação (conhecida e fixa).

    C : list of numpy.ndarray, length M
        Lista de matrizes de autocovariância amostral das inovações.
        C[i] tem shape (nz, nz) e representa Ĉ(i), i = 0, 1, ..., M-1.
        C[0] é a covariância contemporânea (lag 0).

    Retorna
    -------
    Grad : float
        Valor do gradiente da função objetivo .
        
    """
    M = len(C)
    nz, nx = H.shape
    I_nx = np.eye(nx)
    C0 = C[0]

    # Normalização E^2 = diag(C0)^{-1}
    diag_C0 = np.diag(C0)
    diag_C0 = np.maximum(diag_C0, 1e-12)
    E2 = np.diag(1.0 / diag_C0)

    # Matriz de malha fechada
    F_bar = F @ (I_nx - W @ H)

    # Pré‑calcular Phi(i) = H * F_bar^{i-1} * F  para i = 1 .. M-1
    Phi = []
    for i in range(1, M):
        Phi.append(H @ np.linalg.matrix_power(F_bar, i-1) @ F)

    # Calcular X e Z
    X = compute_X(F_bar, F, H, C, M)

    # Inicializar gradiente (nx, nz)
    grad = np.zeros((nx, nz))

    # Somatório principal
    for i in range(1, M):
        # Q_Z = 1/2 * Σ_{j=i}^{M-1} [ Φ(j)^T E^2 C(j) E^2 H + (Φ(j)^T E^2 C(j) E^2 H)^T ]
        Z_i = compute_Z(F_bar, Phi, E2, H, C, start=i, M=M)
        
        # Primeiro termo: Phi(i)^T * E^2 * C(i) * E^2 * C(0)
        term1 = Phi[i-1].T @ E2 @ C[i] @ E2 @ C0
        grad -= term1

        # Segundo termo (soma dupla): l = 0 .. i-2
        for l in range(i-1):
            # [ C(l+1) * E^2 * C(i)^T * E^2 * H * F_bar^{i-l-2} ]^T
            C_l1 = C[l+1]                         # C(ℓ+1)
            power = i - l - 2
            F_pow = np.linalg.matrix_power(F_bar, power)
            inside = C_l1 @ E2 @ C[i].T @ E2 @ H @ F_pow   # (nz, nx)
            term2 = inside.T                      # (nx, nz)
            grad -= term2

        # ---------- Terceiro termo: - F^T Z_i F X (dentro do colchete) ----------
        # - ( - F^T Z_i F X ) = + F^T Z_i F X
        grad -= F.T @ Z_i @ F @ X

    return grad


def gradient_fd(W, F, H, C, eps=1e-5):
    """
    ∂J/∂W calculado por diferenças finitas centradas.
 
    BUG-3 (raiz): a Eq. 60 do artigo deriva ∇J assumindo X = Ψ − W·Ĉ(0)
    (Eq. 57), mas o código estima X via lstsq (Eq. 63). Essas duas
    parametrizações de X são diferentes: o X do lstsq depende implicitamente
    de W através de F̄ e Φ, gerando um gradiente total diferente da Eq. 60.
    O resultado era uma direção de descida incorreta que impedia a convergência.
 
    A solução é calcular ∂J/∂W numericamente, pois o FD captura a dependência
    completa de J em W (inclusive via X recalculado a cada avaliação).
    """
    nx, nz = F.shape[0], H.shape[0]
    grad   = np.zeros((nx, nz))
    for a in range(nx):
        for b in range(nz):
            Wp = W.copy(); Wp[a, b] += eps
            Wm = W.copy(); Wm[a, b] -= eps
            grad[a, b] = (objective_J(Wp, F, H, C) - objective_J(Wm, F, H, C)) / (2 * eps)
    return grad
    
# Método para garantir que F_bar é estável (F de malha fechada)
def is_stable(F_bar):
    """Retorna True se todos autovalores têm magnitude < 1."""
    return np.max(np.abs(np.linalg.eigvals(F_bar))) < 1.0



def update_gain_W(F, H, C, W_init, N, N_s=None, beta=2.0,
                  max_iter=100, tol_W=1e-6, tol_grad=1e-6, tol_J=1e-6,
                  patience=5, c=0.01, c_max=0.2, grad_eps=1e-5):
    """
    Otimiza o ganho W via gradiente descendente com passo adaptativo (bold driver).

    Parâmetros
    ----------
    F : ndarray (nx, nx)
        Matriz de transição de estado.
    H : ndarray (nz, nx)
        Matriz de observação.
    C : list de ndarrays (M, nz, nz)
        Autocovariâncias amostrais das inovações (C[i] = lag i).
    W_init : ndarray (nx, nz)
        Ganho inicial estabilizante.
    N : int
        Número total de amostras observadas (usado para inicializar alpha).
    N_s : int, opcional
        Hiperparâmetro de escala para alpha. Se None, será usado N_s = N.
    beta : float, opcional
        Expoente para ajuste de alpha (padrão 2, conforme artigo).
    max_iter : int
        Máximo de iterações.
    tol : float
        Tolerância para convergência (gradiente e variação de W).
    patience : int
        Iterações sem melhora para parada antecipada.

    Retorna
    -------
    best_W : ndarray (nx, nz)
        Ganho com menor J observado.
    best_J : float
        Valor de J correspondente.
    """
    nx  = F.shape[0]
    I_nx = np.eye(nx)
    if N_s is None:
        N_s = N
 
    alpha  = min(c * (N / N_s) ** beta, c)   # Eq. 136
    c_bar  = min((N / N_s) ** beta, c_max)    # Eq. 138
 
    W       = W_init.copy()
    best_W  = W.copy()
    best_J  = objective_J(W, F, H, C)         # BUG-8: era objective_function
    J_prev  = best_J
    no_improve = 0
 
    for it in range(max_iter):
        grad = gradient_fd(W, F, H, C, eps=grad_eps)   # BUG-3: FD correto
        #grad = gradient_J(W,F,H,C)
        W_new   = W - alpha * grad
        F_bar_n = F @ (I_nx - W_new @ H)
 
        if not is_stable(F_bar_n):
            alpha *= 0.5
            no_improve += 1
            if no_improve >= patience:
                break
            continue
 
        J_new = objective_J(W_new, F, H, C)    # BUG-8
 
        # Bold driver (Eq. 137)
        alpha = 0.5 * alpha if J_new > J_prev else min(1.1 * alpha, c_bar)
 
        W      = W_new
        J_prev = J_new
 
        if J_new < best_J:
            best_J = J_new
            best_W = W.copy()
            no_improve = 0
        else:
            no_improve += 1
 
        # Condições de parada 1-5 (Seção VIII-C1)
        dW = np.linalg.norm(W - W_init) / (np.linalg.norm(W_init) + 1e-12)
        if dW < tol_W:           break   # Cond. 1
        if np.linalg.norm(grad) < tol_grad:  break   # Cond. 2
        if J_new < tol_J:        break   # Cond. 3
        if no_improve >= patience: break # Cond. 4
        # Cond. 5: max_iter já controla o loop
 
    return best_W, best_J


def step3_update_W(v, F, H, W_current, M_lags, max_iter=100, tol=1e-6,
                   N_s=None, beta=2.0):
    """
    Executa o passo 3 completo:
    - Calcula autocovariâncias de v
    - Otimiza W via gradiente adaptativo

    Parâmetros adicionais
    ---------------------
    N_s, beta : parâmetros para inicialização do passo alpha (Eq. 136)
    """
    N = v.shape[0]          # número de amostras observadas
    C = compute_autocovariance(v, M_lags)
    W_opt, J_opt = update_gain_adaptative(F, H, C, W_current, N=N,
                                          N_s=N_s, beta=beta,
                                          max_iter=max_iter, tol=tol)
    return W_opt, J_opt

# Definindo função para estimar R
def estimate_R(nu=None, mu=None, S=None, G=None, H=None, W=None, method='R3'):
    """
    Estima a matriz de covariância de medição R usando um dos cinco métodos do artigo.

    Parâmetros
    ----------
    nu : ndarray (N, nz), opcional
        Sequência de inovações.
    mu : ndarray (N, nz), opcional
        Sequência de resíduos pós-ajuste.
    S : ndarray (nz, nz), opcional
        Covariância amostral das inovações (se já calculada).
    G : ndarray (nz, nz), opcional
        Covariância amostral dos resíduos pós-ajuste (se já calculada).
    H : ndarray (nz, nx)
        Matriz de observação.
    W : ndarray (nx, nz)
        Ganho estacionário do filtro.
    method : str, opcional ('R1', 'R2', 'R3', 'R4', 'R5')
        Método de estimação. Padrão: 'R3' (recomendado).

    Retorna
    -------
    R : ndarray (nz, nz)
        Estimativa da covariância do ruído de medição.
    """
    if S is None and nu is not None:
        S = nu.T @ nu / len(nu)
    if G is None and mu is not None:
        G = mu.T @ mu / len(mu)
 
    I_nz = np.eye(H.shape[0])
 
    if method == 'R1':
        R = (I_nz - H @ W) @ S
    elif method == 'R2':
        if nu is None or mu is None:
            raise ValueError("R2 requer nu e mu.")
        cross = mu.T @ nu / len(nu)
        R = 0.5 * (cross + cross.T)
    elif method == 'R3':
        S_sqrt     = sqrtm(S)
        S_inv_sqrt = inv(S_sqrt)
        R = S_sqrt @ sqrtm(S_inv_sqrt @ G @ S_inv_sqrt) @ S_sqrt
        R = 0.5 * (R + R.T)
    elif method == 'R4':
        R = 0.5 * (G + S - H @ W @ S @ W.T @ H.T)
    elif method == 'R5':
        R = 0.5 * (G @ inv(I_nz - W.T @ H.T) + inv(I_nz - H @ W) @ G)
    else:
        raise ValueError("Método deve ser 'R1'..'R5'.")
 
    R = 0.5 * (R + R.T)
    vals, vecs = np.linalg.eigh(R)
    return np.real(vecs @ np.diag(np.maximum(vals, 1e-12)) @ vecs.T)


# Definição da função do Passo 5
def estimate_Q_and_P(F, H, Gamma, W, R, S, max_iter=100, tol=1e-6, lambda_Q=0.0, mask=None):
    """
    Estima iterativamente Q e P (Passo 5 do artigo).

    Parâmetros
    ----------
    F : ndarray (nx, nx)
        Matriz de transição de estado.
    H : ndarray (nz, nx)
        Matriz de observação.
    Gamma : ndarray (nx, nv)
        Matriz de ganho do ruído de processo.
    W : ndarray (nx, nz)
        Ganho ótimo do filtro (estimado no Passo 3).
    R : ndarray (nz, nz)
        Covariância do ruído de medição (estimada no Passo 4).
    S : ndarray (nz, nz)
        Covariância da inovação (C[0]).
    max_iter : int, opcional
        Máximo de iterações externas para convergência de Q.
    tol : float, opcional
        Tolerância para mudança relativa em Q.
    lambda_Q : float, opcional
        Regularização adicionada a D antes de resolver para Q (Eq. 127).
    mask : ndarray (nv, nv), opcional
        Máscara binária para impor estrutura em Q (ex.: diagonal). 
        Se None, Q é estimado livremente.

    Retorna
    -------
    Q : ndarray (nv, nv)
        Covariância estimada do ruído de processo.
    P : ndarray (nx, nx)
        Covariância de estado atualizada.
    """
    nx   = F.shape[0]
    nv   = Gamma.shape[1]
    I_nx = np.eye(nx)
    F_bar = F @ (I_nx - W @ H)
    GQGt  = W @ S @ W.T        # inicialização Wiener (Eq. 165)
    mask  = np.ones((nv, nv)) if mask is None else mask
 
    Q = P = None
 
    for _ in range(max_iter):
        # Ponto-fixo para P (Eqs. 122-123)
        P_new = solve_discrete_lyapunov(
            F_bar.T,
            W @ R @ W.T + (I_nx - W @ H) @ GQGt @ (I_nx - W @ H).T
        )
        for _ in range(max_iter):
            P_prev = P_new
            P_pred = F @ P_new @ F.T + GQGt
            P_new  = inv(inv(P_pred) + H.T @ inv(R) @ H)
            if np.linalg.norm(P_new - P_prev, 'fro') < tol:
                break
 
        # Atualização de Q (Eqs. 124-127)
        D      = P_new + W @ S @ W.T - F @ P_new @ F.T + lambda_Q * I_nx
        G_pinv = pinv(Gamma)
        Q_new  = G_pinv @ D @ G_pinv.T
        Q_new  = mask * Q_new
        Q_new  = 0.5 * (Q_new + Q_new.T)
        vals, vecs = np.linalg.eigh(Q_new)
        Q_new  = vecs @ np.diag(np.maximum(vals, 1e-12)) @ vecs.T
        GQGt_n = Gamma @ Q_new @ Gamma.T
 
        if Q is not None:
            if np.linalg.norm(Q_new - Q, 'fro') / (np.linalg.norm(Q, 'fro') + 1e-12) < tol:
                Q, P = Q_new, P_new
                break
 
        Q, P, GQGt = Q_new, P_new, GQGt_n
 
    return Q, P

def adaptive_kalman_filter(Z, F, H, Gamma,
                           M_lags=None, max_outer=20, max_inner=100,
                           tol_J=1e-12, tol_W=1e-12, tol_grad=1e-12,
                           lambda_Q=0.0, mask_Q=None,
                           method_R='R3', N_s=None, beta=2.0,
                           patience=5, c=0.01, c_max=0.2,
                           grad_eps=1e-5, verbose=False):
    """
    Executa o algoritmo completo de 6 etapas para estimar W, R, Q e P.

    Parâmetros
    ----------
    Z : ndarray (N, nz)
        Dados de medição.
    F, H, Gamma : ndarray
        Matrizes do sistema linear.
    M_lags : int, opcional
        Número de lags para autocovariância. Se None, usa M = nx.
    max_outer : int
        Máximo de iterações do laço externo (Passo 6).
    max_inner : int
        Máximo de iterações para otimização de W (Passo 3) e Q/P (Passo 5).
    tol_outer : float
        Tolerância para convergência do laço externo.
    tol_inner : float
        Tolerância para otimização interna.
    lambda_Q : float
        Regularização para estimação de Q.
    mask_Q : ndarray (nv, nv), opcional
        Máscara estrutural para Q (ex: np.eye(nv) para diagonal).
    method_R : str
        Método de estimação de R ('R1' a 'R5').
    N_s : int, opcional
        Hiperparâmetro de escala para inicialização do passo alpha (Eq. 136).
        Se None, será usado N_s = N (número de amostras).
    beta : float, opcional
        Expoente para ajuste do passo alpha (padrão 2, conforme artigo).
    verbose : bool
        Se True, imprime progresso a cada iteração.

    Retorna
    -------
    best_W : ndarray (nx, nz)
        Melhor ganho de Kalman encontrado.
    R_est : ndarray (nz, nz)
        Covariância de medição estimada.
    Q_est : ndarray (nv, nv)
        Covariância de processo estimada.
    P_est : ndarray (nx, nx)
        Covariância de estado atualizada.
    history : dict
        Histórico de J, normas de W, R, Q.
    """
    N, nz = Z.shape
    nx    = F.shape[0]
 
    M_lags = max(nx, 5) if M_lags is None else M_lags
 
    W        = init_gain_W(F, H, Gamma)
    best_W   = W.copy()
    best_J   = np.inf
    prev_J   = np.inf
    R_est = Q_est = P_est = None
    history  = {'J': [], 'delta_W': [], 'delta_J': []}
 
    for outer in range(max_outer):
        if verbose:
            print(f"\n--- Iteração externa {outer+1}/{max_outer} ---")
 
        # Passo 1-2: inovações e autocovariâncias com W atual
        nu, _  = generate_residuals(Z, F, H, W)
        C      = compute_autocovariance(nu, M_lags)
 
        # Passo 3: otimizar W
        W_new, J_new = update_gain_W(
            F, H, C, W, N=N, N_s=N_s, beta=beta,
            max_iter=max_inner,
            tol_W=tol_W, tol_grad=tol_grad, tol_J=tol_J,
            patience=patience, c=c, c_max=c_max, grad_eps=grad_eps
        )
 
        if J_new < best_J:
            best_J = J_new
            best_W = W_new.copy()
 
        # BUG-5: regenerar nu/mu com W_new para Steps 4 e 5
        nu_new, mu_new = generate_residuals(Z, F, H, W_new)
        S_new = compute_autocovariance(nu_new, M_lags)[0]
        G_new = mu_new.T @ mu_new / max(N - M_lags, 1)
 
        # Passo 4: estimar R
        R_est = estimate_R(S=S_new, G=G_new, H=H, W=W_new, method=method_R)
 
        # Passo 5: estimar Q e P
        Q_est, P_est = estimate_Q_and_P(
            F, H, Gamma, W_new, R_est, S_new,
            max_iter=max_inner, tol=tol_W,
            lambda_Q=lambda_Q, mask=mask_Q
        )
 
        delta_W = np.linalg.norm(W_new - W, 'fro') / (np.linalg.norm(W, 'fro') + 1e-12)
        delta_J = abs(best_J - prev_J)
        history['J'].append(J_new)
        history['delta_W'].append(delta_W)
        history['delta_J'].append(delta_J)
 
        if verbose:
            print(f"  J={J_new:.6f}  ΔJ={delta_J:.2e}  ΔW={delta_W:.2e}")
 
        W, prev_J = W_new, best_J
 
        if delta_J < tol_J:
            if verbose: print("  Convergência: ΔJ < tol_J")
            break
        if delta_W < tol_W:
            if verbose: print("  Convergência: ΔW < tol_W")
            break
 
    if verbose:
        print(f"\nMelhor J = {best_J:.6f}")
 
    return best_W, R_est, Q_est, P_est, history
 
