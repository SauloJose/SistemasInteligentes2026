# Sistemas Inteligentes — 2026

Visão geral do repositório
- Este workspace reúne materiais, atividades e códigos desenvolvidos para a disciplina "Sistemas Inteligentes" (ano/semestre 2026). O material está organizado por atividades numeradas que cobrem tópicos como filtragem ótima, sintonia de filtros, controle não-linear, controle comportamental, aprendizado por reforço e elaboração de artigos.

Estrutura de pastas (resumo)
- `Artigos/` — Materiais e textos relacionados a artigos ou referências bibliográficas.
- `Atividade 1 - Filtragem Ótima/` — Scripts e modelos MATLAB para introdução à filtragem ótima (veja README nesta pasta).
- `Atividade 2 - Sintonia do Filtro ótimo/` — Pasta para atividades de sintonia (vazia ou com materiais específicos).
- `Atividade 3 - Controle Não-Linear/` — Materiais sobre controle não-linear.
- `Atividade 4  - Controle Comportamental/` — Materiais sobre controle comportamental.
- `Atividade 5 - Reinforced Learning/` — Materiais e experimentos de RL.
- `Atividade 6  - Elaboração de Artigos/` — Apoio para redação e estruturação de artigos.
- `Nova pasta/` — Implementações auxiliares, por exemplo: `EKF.m`, `filtro_kalman_generico.m`, `kalman.m`, `lindberg.m`, `modelsSin.m`.

Dependências e ambiente
- A maior parte dos códigos está em MATLAB (.m). Recomenda-se:
  - MATLAB instalado (ou GNU Octave como alternativa, com possíveis ajustes).
  - Configurar o diretório de trabalho para a pasta raiz do projeto quando for executar scripts que referenciem caminhos relativos.

Como navegar e executar
- Comece por `Atividade 1 - Filtragem Ótima/README.md` para entender os exemplos básicos.
- Abra os arquivos `.m` no MATLAB Editor e execute os scripts diretamente.
- Para reutilizar funções gerais (ex.: `kalman.m` ou `EKF.m`), adicione `Nova pasta/` ao path do MATLAB: `addpath('Nova pasta')`.

Contribuições e organização
- Preserve a estrutura de pastas ao adicionar novos experimentos.
- Nomeie arquivos de forma descritiva (ex.: `atividadeX_experimentoY.m`) e documente o propósito no topo do arquivo.

Contato
- Material preparado para a disciplina; para dúvidas sobre os códigos, consulte o docente ou o responsável pela disciplina.

Notas finais
- Estes READMEs são um ponto de partida — sinta-se à vontade para solicitar que eu gere documentação mais detalhada para qualquer atividade específica (ex.: explicações passo a passo para cada script).