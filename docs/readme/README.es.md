# OpenQuantum Papers

<p align="center"><a href="../../README.md">简体中文</a> · <a href="README.en.md">English</a> · <a href="README.ja.md">日本語</a> · <a href="README.ko.md">한국어</a> · <strong>Español</strong> · <a href="README.fr.md">Français</a> · <a href="README.de.md">Deutsch</a> · <a href="README.pt.md">Português</a> · <a href="README.ru.md">Русский</a> · <a href="README.ar.md">العربية</a></p>

Los seis trabajos de investigación de OpenQuantum abarcan desde una ley de recuperación común para códigos cuánticos aleatorios hasta un simulador de circuitos cuánticos compilado en los pesos de Qwen. Aquí puedes consultar los resultados, las derivaciones y los registros de verificación.

Usa OpenQuantum para recalcular en local el estado fundamental de un sistema de dos cúbits y comparar su energía con una solución exacta independiente. El ejemplo funciona sin clave de modelo.

**[Leer los seis artículos](#papers)** · **[Recalcular el estado fundamental de dos cúbits][openquantum-first-run]** · [Instalar OpenQuantum][openquantum-start]

<a id="papers"></a>

## Seis artículos

Los artículos y sus suplementos están disponibles en PDF en inglés.

### QwenQPU · Un simulador de circuitos cuánticos integrado en Qwen

Cuatro capas estándar de Qwen con pesos fijos ejecutan una simulación clásica por vector de estado de circuitos de seis cúbits. Actualizan las amplitudes complejas y devuelven las probabilidades de medida. Al activar el motor de ejecución en el mismo modelo de lenguaje, las tareas superadas en la evaluación fija pasaron **de 4/48 a 48/48**. El error absoluto máximo de las probabilidades respecto a la referencia de Qiskit fue de **1.42 × 10⁻⁵**.

[Artículo y apéndices][qwenqpu] · [Implementación y datos de evaluación][qwenqpu-code]

![Animación original de QwenQPU. Muestra instrucciones del circuito, 64 direcciones lógicas de amplitudes y las probabilidades de una simulación clásica de seis cúbits.](../../assets/qwenqpu-banner.gif)

*Las amplitudes y la probabilidad final que aparecen a la derecha proceden de una ejecución de las capas estándar de Qwen.*

### Códigos cuánticos aleatorios · Una curva crítica común a la recuperación clásica y cuántica

Para códigos aleatorios obtenidos de un diseño unitario exacto de orden 2, el trabajo establece cotas inferiores y superiores comunes a longitud de bloque finita para los óptimos promediados sobre el conjunto de códigos. Las dos tareas son decodificar el mensaje completo y recuperar el entrelazamiento bajo ruido de Pauli. Para ruido independiente e idénticamente distribuido que cumple las condiciones de regularidad, ambos promedios óptimos se acercan a **la misma curva gaussiana** en la ventana crítica de tasa positiva.

[Artículo][random-codes] · [Suplemento][random-codes-si]

![La banda azul acota los óptimos promediados de ambas tareas. A medida que aumenta la longitud del bloque, la banda se estrecha hacia el límite gaussiano marcado con una línea blanca discontinua.](../../assets/random-code-banner.gif)

*En la figura 1b, para las tres longitudes de bloque, la anchura de las cotas cerca del centro crítico se reduce de unos 31.8 a 3.4 puntos porcentuales.*

### Sistemas disipativos de Majorana · Una correlación a cuatro tiempos selecciona la relajación de la energía

Con interacciones de orden fijo y disipación uniforme, las cuatro inserciones fermiónicas se agrupan en dos pares de intervalo corto. Cada par se suma sobre todas las componentes. Cuando ambos intervalos duran ε, la señal de energía es de orden ε² y los demás modos de decaimiento contribuyen en orden O(ε⁴). El artículo acota el error en una ventana de tiempo finita.

[Artículo][majorana] · [Suplemento][majorana-si]

![La suma de cada par fermiónico, tanto el de excitación como el de lectura, sobre todas las componentes selecciona el modo de energía. A la derecha se muestra la contribución relativa máxima de los otros modos al acortar los intervalos.](../../assets/majorana-banner.gif)

*En los 18 sistemas finitos de la figura 1c, con εJ = 0.16 y gJ = 0.7, la contribución máxima de los otros modos respecto a la señal de energía fue del 0.0983%.*

### Puntos excepcionales de un liouvilliano fermiónico · Dos tasas de decaimiento desde el mismo estado inicial

El número total de partículas y el desequilibrio de ocupación de los modos seleccionan canales de decaimiento distintos desde un mismo estado inicial. Cuando la lectura solo devuelve el sistema retenido o una señal de borrado, los ensayos breves con el estado de vacío alcanzan la ley de escala óptima para el número de preparaciones necesario para distinguir intensidades de bombeo, contando cada preparación.

[Artículo][fermion-ep] · [Suplemento][fermion-ep-si]

![La lectura del número de partículas y del desequilibrio de modos da tasas de decaimiento distintas para el mismo estado inicial. A la derecha se muestran la señal condicionada y la fracción de supervivencia.](../../assets/fermionic-ep-banner.gif)

*En el ejemplo de 12 modos de la figura 1, las tasas de relajación del número total de partículas y del desequilibrio de ocupación son de unos 1.679γ y 1.118γ, respectivamente. Para γt = 5, la fracción de supervivencia es de 1.18 × 10⁻¹⁵.*

### Espectroscopia E₇ en cadenas de Rydberg · Cuatro líneas candidatas con dos sondas

Las sondas uniforme y alternada permiten comparar los canales de paridad par e impar en una misma escala de energía. En una ventana de parámetros del modelo con bloqueo estricto, las razones de masa de cuatro picos candidatos de baja energía se desvían **como máximo un 3%** de las predicciones de E₇. En el punto de trabajo con 36 átomos, la energía de enlace del candidato débilmente ligado A₄ es aproximadamente un 62% mayor que el valor ideal. Esto cuantifica la desviación debida a la longitud finita de la cadena.

[Artículo][e7] · [Suplemento][e7-si]

![Las sondas alternada y uniforme revelan dos picos candidatos de E₇ cada una. A la derecha aparecen los dos espectros y la energía de enlace del estado débilmente ligado para 36 átomos.](../../assets/rydberg-e7-banner.gif)

*En el cálculo con 36 átomos de las figuras 2a y 3a, la energía de enlace del candidato A₄ es el 4.92% de la menor brecha de energía. El valor ideal de E₇ es el 3.04%.*

### Conjuntos abiertos de átomos de Rydberg · Del intercambio de pesos modales a la señal de densidad

Cerca de un punto de pliegue de campo medio, dos modos de decaimiento intercambian pesos en un cruce evitado. Al eliminar las coherencias se recuperan las raíces espectrales y los pesos en el observable. En el ejemplo de 48 átomos con la preparación fija, el peso de un solo modo cambia en unos **8.5 × 10⁻⁴**, y el cambio máximo de la señal de densidad completa es de unos **1.9 × 10⁻⁶**.

[Artículo][rydberg-decay] · [Suplemento][rydberg-decay-si]

![Dos modos de decaimiento intercambian peso espectral con la preparación fija. Una escala logarítmica común compara el cambio de peso de un modo con la diferencia máxima de la señal de densidad completa entre las mismas dos desintonías.](../../assets/rydberg-decay-banner.gif)

*El panel derecho compara el cambio de peso de la figura 5c con la diferencia máxima de la densidad de la figura 5d en las mismas unidades. Entre las dos desintonías marcadas, las magnitudes difieren por un factor de unos 450.*

## Haz tus propios cálculos con OpenQuantum

OpenQuantum es una plataforma de IA de código abierto para investigación cuántica. Está disponible en Codex y como aplicación de escritorio independiente. Puedes simular un circuito cuántico con Qiskit, verificar el resultado por un método independiente y continuar con el mismo código, parámetros y resultados en tu entorno de trabajo.

**[Instalar OpenQuantum][openquantum-start]** · [Ver ejemplos ejecutables][openquantum-examples] · [OpenQuantum en GitHub][openquantum]

## Registros de investigación

[Versiones de los manuscritos y registros de verificación][manuscript-history] · [Archivo de investigación][research-archive] · [Fuentes y figuras del manuscrito de QwenQPU][qwenqpu-tex]

[openquantum]: https://github.com/xi-zhao/OpenQuantum/blob/main/docs/readme/README.es.md
[openquantum-start]: https://github.com/xi-zhao/OpenQuantum/blob/main/docs/readme/README.es.md#quick-start
[openquantum-examples]: https://github.com/xi-zhao/OpenQuantum/blob/main/docs/readme/README.es.md#algorithm-learning-and-method-comparison
[openquantum-first-run]: https://github.com/xi-zhao/OpenQuantum/blob/main/docs/readme/README.es.md#no-model-key-needed-recompute-a-fixed-case
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
