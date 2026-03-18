# Atividade 1 — Filtragem Ótima

Descrição
- Resumo: Esta pasta contém os scripts referentes à Atividade 1 do curso "Sistemas Inteligentes" (2026), focada em filtragem ótima (ex.: filtros de Kalman clássicos e variações). Os exemplos e modelos foram desenvolvidos em MATLAB (.m).

Objetivo da atividade
- Introduzir conceitos de filtragem ótima e implementar modelos/simulações que demonstram o funcionamento do filtro de Kalman e aplicações associadas.

Conteúdo da pasta
- modelo1.m — Script/Modelo 1: exemplo introdutório (configuração básica do sistema, simulação de sinal e aplicação do filtro).
- modelo2.m — Script/Modelo 2: variação do sistema / parâmetros alterados para demonstrar efeitos do ruído.
- modelo3.m — Script/Modelo 3: implementação com medições ruidosas adicionais e análise de convergência.
- modelo4.m — Script/Modelo 4: experimentos com diferentes matrizes de covariância (Q/R) e comparação de desempenho.
- modelo5.m — Script/Modelo 5: complemento com visualizações finais e métricas de erro.

Dependências
- MATLAB (versão compatível com os scripts; estes são escritos em .m puro, sem toolboxes proprietários específicos). Verifique se o caminho de trabalho está apontado para esta pasta.

Como usar
1. Abra o MATLAB e ajuste o diretório de trabalho para esta pasta.
2. Carregue e execute o script desejado, por exemplo:
   - `run('modelo1.m')` ou simplesmente executar o arquivo no Editor/Command Window.
3. Para comparar execuções, altere parâmetros de ruído ou matrizes de covariância dentro dos scripts conforme comentado.

Expectativas de saída
- Gráficos de estado estimado vs. verdadeiro.
- Curvas de erro (erro de estimação, covariâncias) para análise de desempenho.
- Relatórios numéricos simples com RMSE ou outras métricas, dependendo do script.

Boas práticas
- Faça cópias dos arquivos antes de alterar os parâmetros principais.
- Documente alterações locais (comentários ou arquivos separados) quando testar variações experimentais.

Autor / Contato
- Desenvolvido para uso em aulas da disciplina. Consulte o docente/responsável pela disciplina para dúvidas e correções.

Observações
- Há códigos relacionados a EKF e implementações genéricas em outras pastas do projeto; verifique a pasta `Nova pasta` para implementações adicionais (ex.: `EKF.m`).
