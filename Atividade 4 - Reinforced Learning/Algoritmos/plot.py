import pandas as pd
import matplotlib.pyplot as plt
import numpy as np  # Importado para calcular a linha de tendência
import os

# ==========================================
# 1. CARREGAMENTO DOS DADOS
# ==========================================
arquivo_csv = 'Models/dqn/experiment_log.csv'

# Verificação profissional para evitar erros abruptos
if not os.path.exists(arquivo_csv):
    print(f"[ERRO] O arquivo '{arquivo_csv}' não foi encontrado. Verifique o caminho e tente novamente.")
    exit()

df = pd.read_csv(arquivo_csv)

# ==========================================
# 2. INSPEÇÃO DOS DADOS (PRINTS NO TERMINAL)
# ==========================================
pd.set_option('display.max_columns', None)
pd.set_option('display.width', 1000)

print("\n" + "="*60)
print(" 📊 VISÃO GERAL DOS DADOS (PRIMEIRAS LINHAS) ")
print("="*60)
print(df.head().to_string()) 

print("\n" + "="*60)
print(" 📈 RESUMO ESTATÍSTICO ")
print("="*60)
print(df.describe().round(2).to_string()) 
print("="*60 + "\n")

# ==========================================
# 3. CONFIGURAÇÃO E EXIBIÇÃO DOS GRÁFICOS
# ==========================================
print("[INFO] Gerando gráficos. Feche a janela de exibição para encerrar o script...")

plt.figure(figsize=(14, 6))

# --- GRÁFICO 1: Evolução da Recompensa ---
plt.subplot(1, 2, 1)
# Adicionado label='Recompensa' para identificar na legenda
plt.plot(df['episode'], df['reward'], marker='o', color='#1f77b4', linestyle='-', linewidth=1.5, markersize=5, label='Recompensa')

# --- CÁLCULO E PLOT DA LINHA DE TENDÊNCIA ---
# Ajusta uma reta (polinômio de grau 1) aos dados
z = np.polyfit(df['episode'], df['reward'], 1)
p = np.poly1d(z)
# Plota a linha de tendência em vermelho tracejado
plt.plot(df['episode'], p(df['episode']), color='red', linestyle='--', linewidth=2, label='Tendência')

# Customizações do gráfico 1
plt.title('Evolução da Recompensa por Episódio', fontsize=12, fontweight='bold')
plt.xlabel('Episódio')
plt.ylabel('Recompensa (Reward)')
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(loc='upper left') # Ativa a legenda para mostrar 'Recompensa' e 'Tendência'

# --- GRÁFICO 2: Posições do Robô e do Alvo ---
ax2 = plt.subplot(1, 2, 2)

# Plot das posições do robô (Coolwarm: vai do azul [início] para o vermelho [fim])
sc = plt.scatter(df['robot_x'], df['robot_y'], c=df['episode'], cmap='coolwarm', s=50, alpha=0.8, edgecolors='none', label='Robô')

# Pegando as coordenadas iniciais do alvo
target_x = df['target_x'].iloc[0]
target_y = df['target_y'].iloc[0]

# Plot do Alvo (Target)
plt.scatter(target_x, target_y, color='red', marker='X', s=150, zorder=5, label='Alvo (Target)')

# Adicionando o Círculo com Raio 0.5 em volta do alvo
circle = plt.Circle((target_x, target_y), 0.5, color='red', fill=False, linestyle='--', linewidth=1.5, label='Raio de Captura', zorder=4)
ax2.add_patch(circle)

# Força as escalas dos eixos X e Y a serem iguais para o círculo não virar uma elipse
ax2.set_aspect('equal', adjustable='datalim')

# Customizações do gráfico de posições
plt.title('Histórico de Posições no Espaço 2D', fontsize=12, fontweight='bold')
plt.xlabel('Posição X')
plt.ylabel('Posição Y')
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(loc='upper left')

# Barra de cores
cbar = plt.colorbar(sc)
cbar.set_label('Linha do Tempo (Episódios)')

# Exibir
plt.tight_layout()
plt.show()

print("[INFO] Script finalizado com sucesso.")