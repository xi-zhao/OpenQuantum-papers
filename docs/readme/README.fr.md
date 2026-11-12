# OpenQuantum Papers

<p align="center"><a href="../../README.md">简体中文</a> · <a href="README.en.md">English</a> · <a href="README.ja.md">日本語</a> · <a href="README.ko.md">한국어</a> · <a href="README.es.md">Español</a> · <strong>Français</strong> · <a href="README.de.md">Deutsch</a> · <a href="README.pt.md">Português</a> · <a href="README.ru.md">Русский</a> · <a href="README.ar.md">العربية</a></p>

Les six articles de recherche d’OpenQuantum vont d’une loi de récupération commune aux codes quantiques aléatoires à un simulateur de circuits quantiques compilé dans les poids de Qwen. Vous y trouverez les résultats, les démonstrations et les relevés de vérification.

Utilisez OpenQuantum pour recalculer localement l’état fondamental d’un système à deux qubits et comparer son énergie à une solution exacte indépendante. Cet exemple ne requiert aucune clé de modèle.

**[Lire les six articles](#papers)** · **[Recalculer l’état fondamental de deux qubits][openquantum-first-run]** · [Installer OpenQuantum][openquantum-start]

<a id="papers"></a>

## Six articles

Les articles et leurs suppléments sont disponibles en anglais au format PDF.

### QwenQPU · Un simulateur de circuits quantiques intégré à Qwen

Quatre couches standard de Qwen à poids fixes exécutent une simulation classique par vecteur d’état de circuits à six qubits. Elles mettent à jour les amplitudes complexes et calculent les probabilités de mesure. Dans le même modèle de langage, l’activation du moteur d’exécution fait passer le nombre de tâches réussies **de 4/48 à 48/48** sur une évaluation fixe. L’erreur absolue maximale sur les probabilités par rapport à la référence Qiskit est de **1,42 × 10⁻⁵**.

[Article et annexes][qwenqpu] · [Implémentation et données d’évaluation][qwenqpu-code]

![Animation originale de QwenQPU montrant les instructions du circuit, 64 adresses logiques d’amplitudes et les probabilités d’une simulation classique de six qubits.](../../assets/qwenqpu-banner.gif)

*Les amplitudes et la probabilité finale à droite proviennent d’une exécution des couches standard de Qwen.*

### Codes quantiques aléatoires · Une courbe critique commune aux récupérations classique et quantique

Pour des codes aléatoires tirés d’un 2-design unitaire exact, l’article établit des bornes inférieure et supérieure communes à longueur de bloc finie pour la moyenne sur les codes de deux valeurs optimales. Les deux tâches sont le décodage du message entier et la récupération de l’intrication sous bruit de Pauli. Pour un bruit indépendant et identiquement distribué qui satisfait aux conditions de régularité, les deux moyennes tendent vers **la même courbe gaussienne** dans la fenêtre critique à taux de codage positif.

[Article][random-codes] · [Supplément][random-codes-si]

![La bande bleue borne les valeurs optimales moyennées des deux tâches. Quand la longueur de bloc augmente, elle se resserre vers la limite gaussienne représentée par un trait blanc discontinu.](../../assets/random-code-banner.gif)

*Près du centre de la fenêtre critique, la largeur de la bande passe d’environ 31,8 à 3,4 points de pourcentage sur les trois longueurs de bloc de la figure 1b.*

### Systèmes dissipatifs de Majorana · Une corrélation à quatre temps sélectionne la relaxation de l’énergie

Pour des interactions d’ordre fixé et une dissipation uniforme, les quatre insertions fermioniques sont groupées en deux paires à court intervalle. Chaque paire est sommée sur toutes les composantes. Quand les deux intervalles durent ε, le signal d’énergie est d’ordre ε² et la contribution des autres modes de décroissance est en O(ε⁴). L’article borne l’erreur sur une fenêtre temporelle finie.

[Article][majorana] · [Supplément][majorana-si]

![La somme sur toutes les composantes des paires fermioniques de préparation et de lecture sélectionne le mode d’énergie. Le panneau de droite montre la contribution relative maximale des autres modes quand l’intervalle se raccourcit.](../../assets/majorana-banner.gif)

*Pour les 18 systèmes finis de la figure 1c, à εJ = 0,16 et gJ = 0,7, la contribution maximale des autres modes par rapport au signal d’énergie est de 0,0983 %.*

### Points exceptionnels d’un liouvillien fermionique · Deux taux de décroissance pour un même état initial

À partir du même état initial, le nombre total de particules et le déséquilibre d’occupation des modes sélectionnent des canaux de décroissance différents. Lorsque la lecture renvoie le système conservé ou un indicateur d’effacement, des essais courts à partir du vide atteignent la loi d’échelle optimale du nombre de préparations requis pour distinguer des intensités de pompage, en comptant chaque préparation.

[Article][fermion-ep] · [Supplément][fermion-ep-si]

![Le nombre total de particules et le déséquilibre des modes donnent des taux de décroissance distincts pour le même état initial. À droite, le signal conditionnel et la fraction de survie.](../../assets/fermionic-ep-banner.gif)

*Dans l’exemple à 12 modes de la figure 1, les taux de relaxation du nombre total de particules et du déséquilibre valent environ 1,679γ et 1,118γ, respectivement. À γt = 5, la fraction de survie vaut 1,18 × 10⁻¹⁵.*

### Spectroscopie E₇ des chaînes de Rydberg · Quatre raies candidates révélées par deux sondes

Les sondes uniforme et alternée placent les canaux de parité paire et impaire sur une même échelle d’énergie. Dans une fenêtre de paramètres du modèle à blocage strict, les rapports de masse de quatre pics candidats de basse énergie concordent **à 3 % près** avec les prédictions E₇. Au point de fonctionnement à 36 atomes, l’énergie de liaison du candidat faiblement lié A₄ dépasse la valeur idéale d’environ 62 %, ce qui quantifie l’écart dû à la longueur finie de la chaîne.

[Article][e7] · [Supplément][e7-si]

![Les sondes alternée et uniforme révèlent chacune deux pics candidats de E₇. À droite figurent les deux spectres et l’énergie de liaison de l’état faiblement lié pour 36 atomes.](../../assets/rydberg-e7-banner.gif)

*Dans le calcul à 36 atomes des figures 2a et 3a, l’énergie de liaison du candidat A₄ vaut 4,92 % de la plus basse énergie d’excitation. La valeur idéale de E₇ est de 3,04 %.*

### Ensembles ouverts de Rydberg · De l’échange de poids entre modes au signal de densité

Près d’un point de repli du champ moyen, deux modes de décroissance échangent leurs poids à travers un croisement évité. L’élimination des cohérences permet de retrouver les racines spectrales et les poids dans l’observable. Dans l’exemple à 48 atomes à préparation fixée, le poids d’un mode varie d’environ **8,5 × 10⁻⁴**, tandis que la variation maximale du signal de densité complet vaut environ **1,9 × 10⁻⁶**.

[Article][rydberg-decay] · [Supplément][rydberg-decay-si]

![À préparation fixée, deux modes de décroissance échangent du poids spectral. Une échelle logarithmique commune compare la variation du poids d’un mode et l’écart maximal du signal de densité complet entre les deux mêmes désaccords en fréquence.](../../assets/rydberg-decay-banner.gif)

*Le panneau de droite compare la variation du poids du mode de la figure 5c et l’écart maximal de densité de la figure 5d dans la même unité. Entre les deux désaccords indiqués, leurs amplitudes diffèrent d’un facteur d’environ 450.*

## Effectuez vos propres calculs avec OpenQuantum

OpenQuantum est une plateforme d’IA open source pour la recherche quantique, disponible dans Codex et en application de bureau autonome. Simulez un circuit quantique avec Qiskit, vérifiez le résultat indépendamment, puis poursuivez le travail avec le même code, les mêmes paramètres et les résultats dans votre espace de travail.

**[Installer OpenQuantum][openquantum-start]** · [Voir les exemples exécutables][openquantum-examples] · [OpenQuantum sur GitHub][openquantum]

## Dossiers de recherche

[Versions des manuscrits et relevés de vérification][manuscript-history] · [Archives de recherche][research-archive] · [Sources du manuscrit et figures de QwenQPU][qwenqpu-tex]

[openquantum]: https://github.com/xi-zhao/OpenQuantum/blob/main/docs/readme/README.fr.md
[openquantum-start]: https://github.com/xi-zhao/OpenQuantum/blob/main/docs/readme/README.fr.md#quick-start
[openquantum-examples]: https://github.com/xi-zhao/OpenQuantum/blob/main/docs/readme/README.fr.md#algorithm-learning-and-method-comparison
[openquantum-first-run]: https://github.com/xi-zhao/OpenQuantum/blob/main/docs/readme/README.fr.md#no-model-key-needed-recompute-a-fixed-case
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
