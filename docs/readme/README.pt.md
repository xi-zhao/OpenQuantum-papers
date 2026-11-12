# OpenQuantum Papers

<p align="center"><a href="../../README.md">简体中文</a> · <a href="README.en.md">English</a> · <a href="README.ja.md">日本語</a> · <a href="README.ko.md">한국어</a> · <a href="README.es.md">Español</a> · <a href="README.fr.md">Français</a> · <a href="README.de.md">Deutsch</a> · <strong>Português</strong> · <a href="README.ru.md">Русский</a> · <a href="README.ar.md">العربية</a></p>

Os seis artigos de pesquisa do OpenQuantum vão de uma lei de recuperação comum para códigos quânticos aleatórios a um simulador de circuitos quânticos compilado nos pesos do Qwen. Aqui você encontra os resultados, as derivações e os registros de verificação.

Use o OpenQuantum para recalcular localmente o estado fundamental de um sistema de dois qubits e comparar a energia com uma solução exata independente. O exemplo roda sem chave de API do modelo.

**[Ler os seis artigos](#papers)** · **[Recalcular o estado fundamental de dois qubits][openquantum-first-run]** · [Instalar o OpenQuantum][openquantum-start]

<a id="papers"></a>

## Seis artigos

Os artigos e seus suplementos estão disponíveis em PDF em inglês.

### QwenQPU · Um simulador de circuitos quânticos dentro do Qwen

Quatro camadas padrão do Qwen, com pesos fixos, executam uma simulação clássica por vetor de estado de circuitos de seis qubits. Elas atualizam amplitudes complexas e leem as probabilidades de medição. Ao ativar o mecanismo de execução no mesmo modelo de linguagem, o número de tarefas aprovadas na avaliação fixa passou **de 4/48 para 48/48**. O erro absoluto máximo das probabilidades em relação à referência do Qiskit foi de **1,42 × 10⁻⁵**.

[Artigo e apêndices][qwenqpu] · [Implementação e dados da avaliação][qwenqpu-code]

![Animação original do QwenQPU com instruções de circuito, 64 endereços lógicos de amplitudes e probabilidades de uma simulação clássica de seis qubits.](../../assets/qwenqpu-banner.gif)

*As amplitudes e a probabilidade final à direita vêm de uma execução das camadas padrão do Qwen.*

### Códigos quânticos aleatórios · Uma curva crítica comum à recuperação clássica e quântica

Para códigos aleatórios obtidos de um 2-design unitário exato, o artigo estabelece limites inferior e superior comuns em comprimento de bloco finito. Os limites valem para as médias, sobre o conjunto de códigos, dos valores ótimos da decodificação da mensagem inteira e da recuperação do emaranhamento sob ruído de Pauli. Para ruído independente e identicamente distribuído que satisfaz as condições de regularidade, as duas médias se aproximam **da mesma curva gaussiana** na janela crítica de taxa de código positiva.

[Artigo][random-codes] · [Suplemento][random-codes-si]

![A faixa azul delimita os valores ótimos médios das duas tarefas. Com o aumento do comprimento do bloco, ela se estreita em direção ao limite gaussiano indicado pela linha branca tracejada.](../../assets/random-code-banner.gif)

*Perto do centro da janela crítica, a largura da faixa cai de cerca de 31,8 para 3,4 pontos percentuais nos três comprimentos de bloco da figura 1b.*

### Sistemas dissipativos de Majorana · Uma correlação de quatro tempos seleciona a relaxação de energia

Com interações de ordem fixa e dissipação uniforme, as quatro inserções fermiônicas formam dois pares com intervalos curtos. Cada par é somado sobre todas as componentes. Quando os dois intervalos duram ε, o sinal de energia é de ordem ε² e os outros modos de decaimento contribuem em ordem O(ε⁴). O artigo dá um limite para o erro em uma janela de tempo finita.

[Artigo][majorana] · [Suplemento][majorana-si]

![Somar sobre todas as componentes de cada par fermiônico, tanto o de excitação quanto o de leitura, seleciona o modo de energia. À direita está a contribuição relativa máxima dos outros modos à medida que o intervalo diminui.](../../assets/majorana-banner.gif)

*Nos 18 sistemas finitos da figura 1c, com εJ = 0,16 e gJ = 0,7, a contribuição máxima dos outros modos em relação ao sinal de energia foi de 0,0983%.*

### Pontos excepcionais de um liouvilliano fermiônico · Duas taxas de decaimento a partir do mesmo estado inicial

O número total de partículas e o desequilíbrio de ocupação dos modos selecionam canais de decaimento distintos a partir do mesmo estado inicial. Quando a leitura devolve o sistema preservado ou um sinal de apagamento, ensaios curtos com o estado de vácuo atingem a lei de escala ótima do número de preparações necessário para distinguir intensidades de bombeamento, contando cada preparação.

[Artigo][fermion-ep] · [Suplemento][fermion-ep-si]

![O número total de partículas e o desequilíbrio de ocupação apresentam taxas de decaimento distintas a partir do mesmo estado inicial. À direita estão o sinal condicionado e a fração de sobrevivência.](../../assets/fermionic-ep-banner.gif)

*No exemplo de 12 modos da figura 1, as taxas de relaxação do número total de partículas e do desequilíbrio são de aproximadamente 1,679γ e 1,118γ, respectivamente. Em γt = 5, a fração de sobrevivência é de 1,18 × 10⁻¹⁵.*

### Espectroscopia E₇ em cadeias de Rydberg · Quatro linhas candidatas com duas sondas

As sondas uniforme e alternada colocam os canais de paridade par e ímpar na mesma escala de energia. Em uma janela de parâmetros do modelo de bloqueio estrito, as razões de massa de quatro picos candidatos de baixa energia diferem **em até 3%** das previsões de E₇. No ponto de operação com 36 átomos, a energia de ligação do candidato fracamente ligado A₄ é cerca de 62% maior que o valor ideal. Isso quantifica o desvio associado ao comprimento finito da cadeia.

[Artigo][e7] · [Suplemento][e7-si]

![As sondas alternada e uniforme revelam, cada uma, dois picos candidatos de E₇. À direita estão os dois espectros e a energia de ligação do estado fracamente ligado para 36 átomos.](../../assets/rydberg-e7-banner.gif)

*No cálculo com 36 átomos das figuras 2a e 3a, a energia de ligação do candidato A₄ equivale a 4,92% da menor lacuna de energia. O valor ideal de E₇ é 3,04%.*

### Ensembles abertos de Rydberg · Da troca de pesos entre modos ao sinal de densidade

Perto de um ponto de dobra da solução de campo médio, dois modos de decaimento trocam pesos em um cruzamento evitado. Eliminar as coerências permite recuperar as raízes espectrais e os pesos no observável. No exemplo com 48 átomos e preparação fixa, o peso de um único modo varia cerca de **8,5 × 10⁻⁴**, e a maior variação do sinal completo de densidade é de cerca de **1,9 × 10⁻⁶**.

[Artigo][rydberg-decay] · [Suplemento][rydberg-decay-si]

![Com a preparação fixa, dois modos de decaimento trocam peso espectral. Uma escala logarítmica comum compara a variação do peso de um modo com a diferença máxima do sinal completo de densidade entre as mesmas duas dessintonias.](../../assets/rydberg-decay-banner.gif)

*O painel direito compara a variação do peso da figura 5c com a diferença máxima de densidade da figura 5d na mesma unidade. Entre as duas dessintonias marcadas, as magnitudes diferem por um fator de cerca de 450.*

## Faça seus próprios cálculos com o OpenQuantum

OpenQuantum é uma plataforma de IA de código aberto para pesquisa quântica, disponível no Codex e como aplicativo de desktop independente. Simule um circuito quântico com o Qiskit, confira o resultado com um método independente e continue usando o mesmo código, parâmetros e resultados no seu ambiente de trabalho.

**[Instalar o OpenQuantum][openquantum-start]** · [Ver exemplos executáveis][openquantum-examples] · [OpenQuantum no GitHub][openquantum]

## Registros da pesquisa

[Versões dos manuscritos e registros de verificação][manuscript-history] · [Arquivo de pesquisa][research-archive] · [Fontes e figuras do manuscrito do QwenQPU][qwenqpu-tex]

[openquantum]: https://github.com/xi-zhao/OpenQuantum/blob/main/docs/readme/README.pt.md
[openquantum-start]: https://github.com/xi-zhao/OpenQuantum/blob/main/docs/readme/README.pt.md#quick-start
[openquantum-examples]: https://github.com/xi-zhao/OpenQuantum/blob/main/docs/readme/README.pt.md#algorithm-learning-and-method-comparison
[openquantum-first-run]: https://github.com/xi-zhao/OpenQuantum/blob/main/docs/readme/README.pt.md#no-model-key-needed-recompute-a-fixed-case
[qwenqpu]: ../../outputs/qwenqpu/paper/main.pdf
[qwenqpu-code]: https://github.com/xi-zhao/QwenQPU/tree/caae45327b46458163d0840a4c1466e5575a7b4f
[qwenqpu-tex]: https://github.com/xi-zhao/OpenQuantum-papers/tree/1749862c959b44dcc2810319ffeed74a725b08f9/outputs/qwenqpu
[random-codes]: ../../outputs/prl-distinguishability/paper.pdf
[random-codes-si]: ../../outputs/prl-distinguishability/supplement.pdf
[majorana]: ../../outputs/prl-dissipative-syk/paper/main.pdf
[majorana-si]: ../../outputs/prl-dissipative-syk/paper/supplement.pdf
[fermion-ep]: ../../outputs/prl-highorder-ep/paper/manuscript.pdf
[fermion-ep-si]: ../../outputs/prl-highorder-ep/paper/supplement.pdf
[e7]: ../../outputs/e7-rydberg-prl/paper/e7_rydberg.pdf
[e7-si]: ../../outputs/e7-rydberg-prl/paper/e7_rydberg_supplement.pdf
[rydberg-decay]: ../../outputs/prl-chiral-switching/paper/main.pdf
[rydberg-decay-si]: ../../outputs/prl-chiral-switching/paper/supplement.pdf
[manuscript-history]: https://github.com/xi-zhao/OpenQuantum-papers/blob/e9acbf31ee1391e5eb3890450744147d917461c0/README.md
[research-archive]: https://github.com/xi-zhao/OpenQuantum-papers/tree/c/import-research-papers
