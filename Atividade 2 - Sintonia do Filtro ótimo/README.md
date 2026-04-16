# Resultados - Atividade 2: Sintonia do Filtro Ótimo

Este diretório reúne os artefatos usados e gerados durante a Atividade 2, cujo objetivo é implementar, testar e analisar métodos de sintonia de filtros de Kalman adaptativos para sistemas lineares e não lineares.

## Estrutura da pasta

- `sintonia.py` - Implementação central do algoritmo adaptativo de identificação de covariâncias e otimização do ganho de Kalman.
- `utils.py` - Funções utilitárias de suporte, incluindo modelos PV/PVA, geração de trajetórias, filtros EKF e funções matemáticas auxiliares.
- `AdaptAlgorithm.ipynb` - Notebook interativo com análise do algoritmo adaptativo e exemplos de execução.
- `Implementacao.ipynb` - Notebook de implementação que reúne experimentos, simulações e resultados de sintonia.
- `SintoniaKalman.ipynb` - Notebook principal de sintonia do filtro de Kalman usado para testar e visualizar o desempenho do método.
- `Figuras/` - Diretório de figuras geradas automaticamente pelo código para documentação e apresentação de resultados.

## `sintonia.py`

Este arquivo implementa as etapas principais do método de sintonia adaptativa do filtro de Kalman, inspirado no artigo de Zhang et al. (2020). As principais funcionalidades são:

- `minimal_polynomial_coeffs(A, tol=1e-8)`
  - Calcula os coeficientes do polinômio mínimo de `A`.
  - Utiliza uma abordagem de mínimos quadrados para encontrar a relação linear entre potências de `A`.

- `compute_B_G_matrices(F, H, Gamma, W)`
  - Calcula listas de matrizes `B_list` e `G_list` necessárias para construir a matriz de identificabilidade.
  - Estima o polinômio mínimo de `F_bar = F (I - W H)` e gera as matrizes auxiliares conforme o artigo.

- `build_identifiability_matrix(B_list, G_list, m, n_v, n_z)`
  - Constrói a matriz de identificabilidade `I` para verificar se `Q` e `R` são estimáveis a partir dos dados.

- `check_noise_identifiability(F, H, Gamma, W=None, verbose=True)`
  - Avalia se o sistema é identificável para estimar as covariâncias de ruído de processo `Q` e de medição `R`.
  - Retorna indicadores para o caso completo (`Q` simétrica) e para o caso mais restrito (`Q` diagonal).

- `normalized_innovation_squared(nu, S)`
  - Calcula a estatística NIS média para avaliar a consistência das inovações do filtro.

- `init_gain_W(F, H, Gamma, Q0=None, R0=None)`
  - Inicializa um ganho de Kalman `W` estabilizante usando a equação de Riccati discreta.

- `generate_residuals(Z, F, H, W, x0_hat=None)`
  - Executa um filtro de Kalman estacionário sobre a sequência de medições `Z` e gera inovações (`nu`) e resíduos (`mu`).

- `compute_autocovariance(nu, M)`
  - Calcula autocovariâncias amostrais das inovações até `M-1` lags.

- `objective_J(W, F, H, C)`
  - Avalia a função objetivo de sintonia `J(W)` que quantifica correlações residuais da sequência de inovações.

- `lstsq_reg(A, B, reg=1e-6)`
  - Resolve um sistema de mínimos quadrados com regularização de Tikhonov.

- `compute_X(F_bar, F, H, C, M)`
  - Estima a matriz `X` que aparece nas equações de identificação, resolvendo um sistema linear baseado nas autocovariâncias.

- `compute_Z(F_bar, Phi, E2, H, C, M)`
  - Resolve uma equação de Lyapunov discreta para obter a matriz auxiliar `Z` usada no gradiente de `J`.

- `gradient_fd(W, F, H, C, eps=1e-7)`
  - Estima o gradiente de `J` por diferenças finitas centradas.

- `gradient_J(W, F, H, C)`
  - Calcula o gradiente analítico da função objetivo `J` com respeito a `W`.

- `is_stable(F_bar)`
  - Verifica se a matriz de malha fechada `F_bar = F (I - W H)` é estável.

- `update_gain_W(F, H, C, W_init, N, N_s=None, beta=2.0, ...)`
  - Otimiza o ganho de Kalman `W` via gradiente descendente com passo adaptativo (bold driver).
  - Ajusta o passo de aprendizado e garante estabilidade do filtro.

- `step3_update_W(v, F, H, W_current, M_lags, ...)`
  - Executa o passo 3 do algoritmo adaptativo: computa autocovariâncias e atualiza `W`.

- `estimate_R(nu=None, mu=None, S=None, G=None, H=None, W=None, method='R3')`
  - Estima a matriz de covariância de medição `R` usando diferentes métodos (`R1` a `R5`).

- `estimate_Q_and_P(F, H, Gamma, W, R, S, ...)`
  - Estima iterativamente `Q` e `P` conforme as equações do algoritmo adaptativo.
  - Aplica regularização e garante que `Q` seja simétrica definida positiva.

- `adaptive_kalman_filter(Z, F, H, Gamma, ...)`
  - Implementa o ciclo completo do algoritmo adaptativo de sintonia:
    1. gera inovações
    2. calcula autocovariâncias
    3. otimiza `W`
    4. estima `R`
    5. estima `Q` e `P`
    6. repete até convergência.

- `run_adaptive_kalman_experiment(F, H, Gamma, Q_true, R_true, ...)`
  - Executa um experimento completo com geração sintética de dados.
  - Gera as trajetórias, executa o algoritmo adaptativo e produz gráficos de convergência e análises NIS.

## `utils.py`

O módulo `utils.py` traz funções de suporte para modelos de movimento e filtragem, com foco em localização baseada em distância e filtros estendidos.

### Funções utilitárias gerais

- `dist(point_a, point_b)`
  - Calcula distância euclidiana entre dois pontos.

- `h_obs(pos, bases)`
  - Calcula o vetor de observação `h(x)` como distâncias entre a posição `pos` e as estações base.

- `create_cov(sigmas)`
  - Constrói uma matriz de covariância diagonal a partir de um vetor de desvios padrão.

- `get_noise_vec(R)`
  - Gera um vetor de ruído branco com covariância `R`.

### Modelo PV (Posição-Velocidade)

- `get_Q_disc_PV(q_diag, dt)`
  - Gera a matriz de covariância de processo discreta `Q` para o modelo PV.

- `calc_h_and_H_PV(x, bases)`
  - Calcula a medição esperada `h(x)` e a matriz Jacobiana `H` para o modelo PV.

- `EKF_PV`
  - Classe que implementa o filtro de Kalman estendido para o modelo PV.
  - Métodos:
    - `predict()` - prediz estado e covariância.
    - `update(z, bases)` - corrige o estado usando a medição `z`.

### Modelo PVA (Posição-Velocidade-Aceleração)

- `get_Q_disc_PVA(q_diag, dt)`
  - Calcula a matriz de covariância de processo discreta `Q` para o modelo PVA.

- `calc_h_and_H_PVA(x, bases)`
  - Calcula a medição esperada `h(x)` e Jacobiana `H` para o modelo PVA.

- `EKF_PVA`
  - Classe que implementa o filtro de Kalman estendido para o modelo PVA.
  - Métodos:
    - `predict()` - prediz o próximo estado.
    - `update(z, bases)` - corrige a estimativa com a medição.

### Geração de trajetórias e visualização

- `TrajectoryGenerator`
  - Classe para gerar trajetórias de referência e visualizar cenários.
  - Métodos:
    - `_init_vectors(num_steps)` - inicializa vetores de tempo e estado.
    - `generate_circle(radius, center, linear_velocity, duration)` - gera uma trajetória circular com cinemática completa.
    - `generate_square(side_length, bottom_left, linear_velocity)` - gera uma trajetória quadrada com segmentos lineares.
    - `generate_tanh_curve(start_pos, end_pos, amplitude_y, smoothness, duration)` - gera uma curva suave em S usando função tangente hiperbólica.
    - `plot_scenario(t, states_gt, title)` - plota a trajetória e as variáveis cinemáticas.

## Notebooks

- `AdaptAlgorithm.ipynb`
  - Notebook destinado à análise do algoritmo adaptativo.
  - Geralmente utilizado para testar parâmetros, visualizar convergência e inspecionar resultados intermediários.

- `Implementacao.ipynb`
  - Notebook de implementação prática.
  - Normalmente contém exemplos de execução do filtro, experimentos e comparações de desempenho.

- `SintoniaKalman.ipynb`
  - Notebook principal para a atividade de sintonia do filtro de Kalman.
  - Usado para rodar experimentos interativos, ajustar parâmetros de covariância e gerar gráficos explicativos.

## Recomendações de uso

1. Execute `sintonia.py` ou os notebooks para carregar o algoritmo adaptativo.
2. Utilize `utils.py` para simular trajetórias e gerar medições ruidosas.
3. Avalie a identificabilidade com `check_noise_identifiability(...)` antes de estimar `Q` e `R`.
4. Ajuste o número de lags `M_lags`, a regularização `lambda_Q` e o método `method_R` conforme o problema.

## Observações

- O código contém funções voltadas para pesquisa e análise numérica de filtros adaptativos.
- O arquivo `sintonia.py` inclui implementação completa do algoritmo de estimação iterativa de `W`, `R`, `Q` e `P`.
- O `README` é um ponto de partida para compreender a estrutura do diretório e a finalidade de cada função/documento.
