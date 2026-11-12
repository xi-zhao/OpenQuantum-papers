# OpenQuantum Papers

<p align="center"><a href="../../README.md">简体中文</a> · <a href="README.en.md">English</a> · <a href="README.ja.md">日本語</a> · <strong>한국어</strong> · <a href="README.es.md">Español</a> · <a href="README.fr.md">Français</a> · <a href="README.de.md">Deutsch</a> · <a href="README.pt.md">Português</a> · <a href="README.ru.md">Русский</a> · <a href="README.ar.md">العربية</a></p>

OpenQuantum이 완성한 연구 논문 여섯 편입니다. 무작위 양자 코드에 공통으로 성립하는 복구 법칙부터 Qwen의 가중치에 내장한 양자 회로 시뮬레이터까지, 결과와 유도 과정, 검증 기록을 확인할 수 있습니다.

OpenQuantum으로 2큐비트 계의 바닥상태를 로컬에서 다시 계산하고, 에너지를 독립적인 정확해와 비교해 보세요. 모델 API 키 없이 실행할 수 있습니다.

**[논문 여섯 편 읽기](#papers)** · **[2큐비트 바닥상태 다시 계산하기][openquantum-first-run]** · [OpenQuantum 설치][openquantum-start]

<a id="papers"></a>

## 논문 여섯 편

논문과 보충 자료는 영어 PDF로 제공됩니다.

### QwenQPU · Qwen에 내장한 양자 회로 시뮬레이터

가중치를 고정한 표준 Qwen 레이어 네 개가 6큐비트 회로를 고전적인 상태 벡터 방식으로 시뮬레이션합니다. 복소 진폭을 갱신하고 측정 확률을 읽어 냅니다. 같은 언어 모델에서 실행기를 켜자 고정된 평가 과제 48개 중 통과 건수가 **4/48에서 48/48로** 늘었습니다. Qiskit 참조값에 대한 확률의 최대 절대 오차는 **1.42 × 10⁻⁵**였습니다.

[논문 및 부록][qwenqpu] · [구현 및 평가 데이터][qwenqpu-code]

![QwenQPU 원본 배너. 회로 명령, 논리 진폭 주소 64개, 6큐비트 고전 시뮬레이션의 확률 출력을 보여 줍니다.](../../assets/qwenqpu-banner.gif)

*오른쪽의 진폭과 최종 확률은 표준 Qwen 레이어의 실행 기록에서 가져왔습니다.*

### 무작위 양자 코드 · 고전 정보와 양자 정보 복구의 공통 임계 곡선

정확한 유니터리 2-디자인에서 생성한 무작위 코드에 대해, 파울리 잡음 아래의 전체 메시지 복호화와 얽힘 복구를 다룹니다. 두 최적값의 앙상블 평균에 공통으로 적용되는 유한 블록 길이의 상한과 하한을 얻었습니다. 정칙성 조건을 만족하는 독립 동일 분포 잡음에서는 양의 부호율 임계 구간에서 두 평균이 **같은 가우스 곡선**에 가까워집니다.

[논문][random-codes] · [보충 자료][random-codes-si]

![파란 띠가 두 과제의 최적값을 앙상블 평균한 값을 모두 감쌉니다. 블록 길이가 늘어나면 흰색 점선으로 표시한 가우스 극한을 향해 띠가 좁아집니다.](../../assets/random-code-banner.gif)

*논문 그림 1b에서 세 가지 블록 길이를 비교하면, 임계 구간 중심 부근의 경계 폭이 약 31.8%p에서 3.4%p로 줄어듭니다.*

### 소산 마요라나 계 · 네 시각 상관함수로 에너지 완화 선택하기

고정 차수의 상호작용과 균일한 소산 아래에서 네 번의 페르미온 삽입을 짧은 시간 간격의 두 쌍으로 묶고, 각 쌍을 모든 성분에 대해 합합니다. 두 쌍의 시간 간격이 모두 ε일 때 에너지 신호는 ε² 차수이며, 다른 감쇠 모드의 기여는 O(ε⁴)입니다. 논문은 유한 시간 구간에서 오차의 상한을 제시합니다.

[논문][majorana] · [보충 자료][majorana-si]

![여기와 판독에 쓰이는 페르미온 쌍을 각각 모든 성분에 대해 합해 에너지 모드를 선택합니다. 오른쪽은 쌍의 시간 간격을 줄일 때 다른 모드가 기여하는 최대 상대 비율을 보여 줍니다.](../../assets/majorana-banner.gif)

*그림 1c의 유한계 18개에서 εJ = 0.16, gJ = 0.7일 때 에너지 신호에 대한 다른 모드의 상대 기여는 최대 0.0983%였습니다.*

### 페르미온 리우빌리안 예외점 · 같은 초기 상태에서 나타나는 두 감쇠율

같은 초기 상태라도 총 입자 수와 모드 점유수 불균형은 서로 다른 감쇠 채널을 선택합니다. 판독 결과가 유지된 계 또는 소거 표지로 한정될 때 모든 상태 준비를 비용에 포함하면, 진공 상태를 사용한 짧은 시간의 시행이 펌핑 세기를 구별하는 데 필요한 최적의 준비 횟수 스케일링을 달성합니다.

[논문][fermion-ep] · [보충 자료][fermion-ep-si]

![같은 초기 상태에서 총 입자 수와 모드 점유수 불균형을 읽으면 서로 다른 감쇠율이 나옵니다. 오른쪽은 조건부 신호와 생존 비율입니다.](../../assets/fermionic-ep-banner.gif)

*그림 1의 12모드 사례에서 총 입자 수와 모드 점유수 불균형의 완화율은 각각 약 1.679γ와 1.118γ입니다. γt = 5에서 생존 비율은 1.18 × 10⁻¹⁵입니다.*

### 리드베리 사슬의 E₇ 분광 · 두 탐침에서 포착한 네 후보선

균일 탐침과 교대 탐침으로 짝수 및 홀수 패리티 채널을 같은 에너지 척도에 놓았습니다. 봉쇄 조건을 엄밀히 적용한 모델의 한 매개변수 구간에서 저에너지 후보 피크 네 개의 질량비는 E₇ 예측값과 **3% 이내**로 일치합니다. 원자 36개인 계산점에서 얕은 결합 상태 A₄ 후보의 결합 에너지는 이상값보다 약 62% 높아, 유한 사슬의 편차를 정량화할 수 있습니다.

[논문][e7] · [보충 자료][e7-si]

![교대 탐침과 균일 탐침으로 각각 E₇ 후보 피크 두 개를 확인합니다. 오른쪽은 원자 36개에서 얻은 두 채널의 스펙트럼과 얕은 결합 상태의 결합 에너지입니다.](../../assets/rydberg-e7-banner.gif)

*그림 2a와 3a의 원자 36개 계산에서 A₄ 후보의 결합 에너지는 가장 낮은 에너지 갭의 4.92%입니다. E₇의 이상값은 3.04%입니다.*

### 열린 리드베리 앙상블 · 감쇠 모드의 가중치 교환에서 밀도 신호까지

평균장의 접힘점 근처에서 두 감쇠 모드가 회피 교차를 거치며 가중치를 교환합니다. 결맞음 성분을 소거해 스펙트럼 방정식의 근과 관측량에 대한 가중치를 복원합니다. 준비 조건을 고정한 원자 48개의 사례에서 개별 모드의 가중치는 약 **8.5 × 10⁻⁴** 변했으며, 전체 밀도 신호의 최대 변화는 약 **1.9 × 10⁻⁶**였습니다.

[논문][rydberg-decay] · [보충 자료][rydberg-decay-si]

![준비 조건을 고정했을 때 두 감쇠 모드가 스펙트럼 가중치를 교환합니다. 같은 로그 척도에서 두 디튜닝 값 사이의 개별 모드 가중치 변화와 전체 밀도 신호의 최대 차이를 비교합니다.](../../assets/rydberg-decay-banner.gif)

*오른쪽은 그림 5c의 모드 가중치 변화와 그림 5d의 전체 밀도 최대 차이를 같은 밀도 단위로 나타냅니다. 표시된 두 디튜닝 값에서 두 변화의 크기는 약 450배 차이가 납니다.*

## OpenQuantum으로 직접 계산하기

OpenQuantum은 Codex 또는 독립 데스크톱 앱에서 사용할 수 있는 양자 연구용 오픈소스 AI 플랫폼입니다. 양자 회로를 Qiskit으로 시뮬레이션하고 독립적인 방법으로 결과를 확인한 뒤, 같은 작업 공간의 코드, 매개변수, 출력을 이어서 활용할 수 있습니다.

**[OpenQuantum 설치][openquantum-start]** · [실행 가능한 계산 예제 보기][openquantum-examples] · [GitHub의 OpenQuantum][openquantum]

## 연구 기록

[원고 버전 및 검증 기록][manuscript-history] · [연구 아카이브][research-archive] · [QwenQPU 원고 소스 및 그림][qwenqpu-tex]

[openquantum]: https://github.com/xi-zhao/OpenQuantum/blob/main/docs/readme/README.ko.md
[openquantum-start]: https://github.com/xi-zhao/OpenQuantum/blob/main/docs/readme/README.ko.md#quick-start
[openquantum-examples]: https://github.com/xi-zhao/OpenQuantum/blob/main/docs/readme/README.ko.md#algorithm-learning-and-method-comparison
[openquantum-first-run]: https://github.com/xi-zhao/OpenQuantum/blob/main/docs/readme/README.ko.md#no-model-key-needed-recompute-a-fixed-case
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
