import matplotlib.pyplot as plt

def plot_rapido(x, y=None, titulo="Meu Gráfico", xlabel="Eixo X", ylabel="Eixo Y", tamanho=(10, 6), grid=True, **kwargs):
    """
    Wrapper rápido para gerar gráficos 2D padronizados no Matplotlib.
    """
    # 1. Configura o tamanho
    plt.figure(figsize=tamanho)
    
    # 2. Plota os dados (aceita só Y, ou X e Y)
    if y is None:
        plt.plot(x, **kwargs)
    else:
        plt.plot(x, y, **kwargs)
        
    # 3. Textos e Labels (com um pequeno ajuste de fonte para ficar mais elegante)
    plt.title(titulo, fontsize=14, fontweight='bold')
    plt.xlabel(xlabel, fontsize=12)
    plt.ylabel(ylabel, fontsize=12)
    
    # 4. Grid opcional estilizado
    if grid:
        plt.grid(True, linestyle='--', alpha=0.6)
        
    # Ajusta as margens e exibe
    plt.tight_layout()
    plt.show()