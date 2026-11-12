# OpenQuantum Papers

<p align="center"><a href="../../README.md">简体中文</a> · <a href="README.en.md">English</a> · <a href="README.ja.md">日本語</a> · <a href="README.ko.md">한국어</a> · <a href="README.es.md">Español</a> · <a href="README.fr.md">Français</a> · <strong>Deutsch</strong> · <a href="README.pt.md">Português</a> · <a href="README.ru.md">Русский</a> · <a href="README.ar.md">العربية</a></p>

Die sechs Forschungsarbeiten von OpenQuantum reichen von einem gemeinsamen Wiederherstellungsgesetz für zufällige Quantencodes bis zu einem Quantenschaltkreis-Simulator, der in die Gewichte von Qwen einkompiliert ist. Hier finden Sie die Ergebnisse, Herleitungen und Prüfprotokolle.

Berechnen Sie mit OpenQuantum den Grundzustand eines Zwei-Qubit-Systems lokal nach und vergleichen Sie die Energie mit einer unabhängigen exakten Lösung. Für das Beispiel ist kein API-Schlüssel für ein Modell erforderlich.

**[Die sechs Arbeiten lesen](#papers)** · **[Den Zwei-Qubit-Grundzustand nachrechnen][openquantum-first-run]** · [OpenQuantum installieren][openquantum-start]

<a id="papers"></a>

## Sechs Forschungsarbeiten

Die Arbeiten und ihre Ergänzungsmaterialien liegen als englische PDFs vor.

### QwenQPU · Ein Quantenschaltkreis-Simulator in Qwen

Vier Standardschichten von Qwen mit festen Gewichten simulieren Schaltkreise mit sechs Qubits klassisch über ihren Zustandsvektor. Sie aktualisieren komplexe Amplituden und lesen Messwahrscheinlichkeiten aus. Mit aktivierter Schaltkreisausführung stieg die Zahl der bestandenen Aufgaben im selben Sprachmodell **von 4/48 auf 48/48**. Der maximale absolute Fehler der Wahrscheinlichkeiten gegenüber der Qiskit-Referenz betrug **1,42 × 10⁻⁵**.

[Arbeit und Anhänge][qwenqpu] · [Implementierung und Evaluationsdaten][qwenqpu-code]

![Originalanimation von QwenQPU mit Schaltkreisbefehlen, 64 logischen Amplitudenadressen und Wahrscheinlichkeiten aus einer klassischen Simulation mit sechs Qubits.](../../assets/qwenqpu-banner.gif)

*Die rechts dargestellten Amplituden und die abschließende Wahrscheinlichkeit stammen aus einem Lauf der Standardschichten von Qwen.*

### Zufällige Quantencodes · Eine gemeinsame kritische Kurve für klassische und Quanteninformation

Für zufällige Codes aus einem exakten unitären 2-Design gibt die Arbeit gemeinsame obere und untere Schranken bei endlicher Blocklänge an. Die Schranken gelten für die über das Code-Ensemble gemittelten Optima der Decodierung der vollständigen Nachricht und der Wiederherstellung von Verschränkung unter Pauli-Rauschen. Bei unabhängig und identisch verteiltem Rauschen unter den angegebenen Regularitätsbedingungen nähern sich beide Mittelwerte im kritischen Fenster bei positiver Coderate **derselben Gauß-Kurve**.

[Arbeit][random-codes] · [Ergänzungsmaterial][random-codes-si]

![Das blaue Band begrenzt die Ensemble-Mittelwerte der Optima beider Aufgaben. Mit wachsender Blocklänge verengt es sich zur weiß gestrichelten Gauß-Grenzkurve hin.](../../assets/random-code-banner.gif)

*In Abb. 1b sinkt die Breite des Bands nahe der Mitte des kritischen Fensters über drei Blocklängen hinweg von etwa 31,8 auf 3,4 Prozentpunkte.*

### Dissipative Majorana-Systeme · Eine Vierzeitkorrelation wählt die Energierelaxation aus

Bei Wechselwirkungen fester Ordnung und gleichförmiger Dissipation bilden die vier Fermioneneinfügungen zwei Paare mit kurzen Zeitabständen. Jedes Paar wird über alle Komponenten summiert. Haben beide Paare den Abstand ε, so hat das Energiesignal die Ordnung ε² und der Beitrag der übrigen Zerfallsmoden die Ordnung O(ε⁴). Die Arbeit gibt eine Fehlerschranke für ein endliches Zeitfenster an.

[Arbeit][majorana] · [Ergänzungsmaterial][majorana-si]

![Die Summation der Fermionenpaare bei Anregung und Auslese über jeweils alle Komponenten wählt die Energiemode aus. Rechts ist der maximale relative Beitrag der übrigen Moden für kürzere Paarabstände zu sehen.](../../assets/majorana-banner.gif)

*In den 18 endlichen Systemen aus Abb. 1c betrug bei εJ = 0,16 und gJ = 0,7 der maximale Beitrag der übrigen Moden relativ zum Energiesignal 0,0983 %.*

### Exzeptionelle Punkte eines fermionischen Liouvillians · Zwei Zerfallsraten aus demselben Anfangszustand

Gesamtteilchenzahl und Besetzungsungleichgewicht der Moden wählen aus demselben Anfangszustand unterschiedliche Zerfallskanäle aus. Wenn die Auslese das erhaltene System oder ein Erasure-Flag liefert, erreichen kurze Vakuumversuche die optimale Skalierung der nötigen Zahl an Präparationen zur Unterscheidung von Pumpstärken. Dabei zählt jede Präparation.

[Arbeit][fermion-ep] · [Ergänzungsmaterial][fermion-ep-si]

![Aus demselben Anfangszustand liefern Teilchenzahl und Besetzungsungleichgewicht verschiedene Zerfallsraten. Rechts stehen das bedingte Signal und der Überlebensanteil.](../../assets/fermionic-ep-banner.gif)

*Im Beispiel mit 12 Moden in Abb. 1 betragen die Relaxationsraten von Gesamtteilchenzahl und Besetzungsungleichgewicht etwa 1,679γ beziehungsweise 1,118γ. Bei γt = 5 beträgt der Überlebensanteil 1,18 × 10⁻¹⁵.*

### E₇-Spektroskopie in Rydberg-Ketten · Vier Linienkandidaten aus zwei Sonden

Homogene und alternierende Sonden bringen die Kanäle gerader und ungerader Parität auf dieselbe Energieskala. In einem Parameterfenster des strikten Blockademodells weichen die Massenverhältnisse von vier niederenergetischen Peak-Kandidaten um **höchstens 3 %** von den E₇-Vorhersagen ab. Am Arbeitspunkt mit 36 Atomen liegt die Bindungsenergie des schwach gebundenen Kandidaten A₄ etwa 62 % über dem Idealwert. Das quantifiziert die Abweichung für die endliche Kette.

[Arbeit][e7] · [Ergänzungsmaterial][e7-si]

![Alternierende und homogene Sonden zeigen jeweils zwei E₇-Peak-Kandidaten. Rechts sind die beiden Spektren und die Bindungsenergie des schwach gebundenen Zustands für 36 Atome zu sehen.](../../assets/rydberg-e7-banner.gif)

*In der Rechnung mit 36 Atomen in Abb. 2a und 3a beträgt die Bindungsenergie des A₄-Kandidaten 4,92 % der niedrigsten Energielücke. Der ideale E₇-Wert beträgt 3,04 %.*

### Offene Rydberg-Ensembles · Vom Gewichtsaustausch der Zerfallsmoden zum Dichtesignal

In der Nähe eines Faltpunkts der Mittelfeldlösung tauschen zwei Zerfallsmoden an einer vermiedenen Kreuzung ihre Gewichte. Durch Elimination der Kohärenzen lassen sich die Nullstellen der Spektralgleichung und die Gewichte in der Observablen rekonstruieren. Im Beispiel mit 48 Atomen bei fester Präparation ändert sich das Gewicht einer Mode um etwa **8,5 × 10⁻⁴**. Die größte Änderung des vollständigen Dichtesignals beträgt etwa **1,9 × 10⁻⁶**.

[Arbeit][rydberg-decay] · [Ergänzungsmaterial][rydberg-decay-si]

![Zwei Zerfallsmoden tauschen bei fester Präparation spektrales Gewicht aus. Auf einer gemeinsamen logarithmischen Skala stehen die Gewichtsänderung einer Mode und die größte Differenz des vollständigen Dichtesignals zwischen denselben zwei Verstimmungen.](../../assets/rydberg-decay-banner.gif)

*Rechts werden die Gewichtsänderung aus Abb. 5c und die größte Dichtedifferenz aus Abb. 5d in derselben Dichteeinheit verglichen. Zwischen den markierten Verstimmungen unterscheiden sich die Beträge um etwa den Faktor 450.*

## Eigene Berechnungen mit OpenQuantum

OpenQuantum ist eine Open-Source-KI-Plattform für die Quantenforschung und läuft in Codex oder als eigenständige Desktop-App. Simulieren Sie einen Quantenschaltkreis mit Qiskit, prüfen Sie das Ergebnis unabhängig und arbeiten Sie in derselben Workbench mit dem Code, den Parametern und den Ausgaben weiter.

**[OpenQuantum installieren][openquantum-start]** · [Ausführbare Beispiele ansehen][openquantum-examples] · [OpenQuantum auf GitHub][openquantum]

## Forschungsaufzeichnungen

[Manuskriptversionen und Prüfprotokolle][manuscript-history] · [Forschungsarchiv][research-archive] · [QwenQPU-Manuskriptquellen und Abbildungen][qwenqpu-tex]

[openquantum]: https://github.com/xi-zhao/OpenQuantum/blob/main/docs/readme/README.de.md
[openquantum-start]: https://github.com/xi-zhao/OpenQuantum/blob/main/docs/readme/README.de.md#quick-start
[openquantum-examples]: https://github.com/xi-zhao/OpenQuantum/blob/main/docs/readme/README.de.md#algorithm-learning-and-method-comparison
[openquantum-first-run]: https://github.com/xi-zhao/OpenQuantum/blob/main/docs/readme/README.de.md#no-model-key-needed-recompute-a-fixed-case
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
