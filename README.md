# Projeto: Sistemas Inteligentes 2026

Este repositório contém o desenvolvimento de atividades relacionadas à filtragem ótima, controle não-linear e diagnóstico de sistemas. O objetivo principal é implementar e analisar algoritmos avançados para controle e filtragem em sistemas dinâmicos, utilizando ferramentas como MATLAB e Python.

---

## Funcionalidades do Projeto

O projeto está dividido em três principais atividades, cada uma com objetivos específicos:

1. **Filtragem Ótima**: Implementação de filtros de Kalman e suas variantes para estimativa de estados em sistemas dinâmicos.
2. **Sintonia de Filtros Ótimos**: Ajuste de parâmetros e análise de desempenho de filtros ótimos.
3. **Controle Não-Linear**: Desenvolvimento de estratégias de controle e diagnóstico para sistemas não-lineares.

Cada atividade está organizada em pastas específicas, contendo scripts, modelos e resultados.

---

## Estrutura do Repositório

Abaixo está a descrição detalhada da estrutura do repositório, incluindo as pastas e os arquivos principais:

### **1. Atividade 1 - Filtragem Ótima**
- **Descrição**: Contém a implementação de filtros de Kalman e suas variantes (ex.: Filtro de Kalman Estendido - EKF).
- **Arquivos**:
  - `modelo1.m`: Script MATLAB para simulação de um sistema linear com filtro de Kalman.
  - `modelo2.m`: Script MATLAB para simulação de um sistema não-linear com EKF.
  - `kalman.m`: Implementação genérica do filtro de Kalman.
  - `RESULTADOS/`: Pasta contendo os resultados das simulações, como gráficos e dados gerados.
- **Utilização**:
  - Execute os scripts no MATLAB para simular os sistemas e gerar os resultados.
  - Analise os gráficos gerados para avaliar o desempenho dos filtros.

---

### **2. Atividade 2 - Sintonia do Filtro Ótimo**
- **Descrição**: Foco na sintonia de parâmetros dos filtros ótimos para melhorar o desempenho.
- **Arquivos**:
  - `sintonia_kalman.m`: Script MATLAB para ajuste de parâmetros do filtro de Kalman.
  - `analise_resultados.m`: Script para análise dos resultados obtidos após a sintonia.
  - `Artigos/`: Referências teóricas utilizadas como base para a implementação.
  - `RESULTADOS/`: Dados e gráficos gerados durante os experimentos.
- **Utilização**:
  - Ajuste os parâmetros no script `sintonia_kalman.m` e execute no MATLAB.
  - Utilize `analise_resultados.m` para visualizar o impacto das alterações.

---

### **3. Atividade 3 - Controle Não-Linear**
- **Descrição**: Desenvolvimento de estratégias de controle e diagnóstico para sistemas não-lineares.
- **Arquivos**:
  - `Diagnostico.ipynb`: Notebook Python com análises e simulações de controle não-linear.
  - `utils.py`: Funções auxiliares para cálculos e simulações.
  - `TE_process/`: Scripts relacionados ao processo TE (Tennessee Eastman).
- **Utilização**:
  - Abra o notebook `Diagnostico.ipynb` no Jupyter Notebook e execute as células para realizar as análises.
  - Utilize as funções em `utils.py` para suporte às simulações.

---

### **4. Nova Pasta**
- **Descrição**: Scripts adicionais para implementação de filtros e modelos.
- **Arquivos**:
  - `EKF.m`: Implementação do Filtro de Kalman Estendido.
  - `filtro_kalman_generico.m`: Script genérico para filtros de Kalman.
- **Utilização**:
  - Scripts prontos para serem adaptados e utilizados em novos projetos.

---

## Requisitos do Projeto

Para executar os scripts e notebooks deste projeto, é necessário ter as seguintes ferramentas instaladas:

### **MATLAB**
- Versão recomendada: R2021a ou superior.
- Necessário para executar os scripts `.m` das atividades 1 e 2.

### **Python**
- Versão recomendada: 3.8 ou superior.
- Bibliotecas necessárias:
  - `numpy`
  - `matplotlib`
  - `pandas`
  - `scipy`
- Utilize o comando abaixo para instalar as dependências:
  ```bash
  pip install numpy matplotlib pandas scipy
  ```
---

## Como Executar
1. Filtragem Ótima e Sintonia:
- Abra o MATLAB e navegue até a pasta correspondente.
- Execute os scripts .m para realizar as simulações.
2. Controle Não-Linear:
- Abra o arquivo `Diagnostico.ipynb` no Jupyter Notebook.
- Execute as células sequencialmente para realizar as análises.
### Contribuição
Contribuições são bem-vindas! Para contribuir:

1. Faça um fork deste repositório.
2. Crie uma branch para sua funcionalidade: `git checkout -b minha-funcionalidade`.
3. Envie um pull request com suas alterações.

### Licença
Este projeto está licenciado sob a licença MIT. Consulte o arquivo LICENSE para mais detalhes.

### Contato
Para dúvidas ou sugestões, entre em contato pelo e-mail: `seuemail@exemplo.com.`