# OpenQuantum Papers

<p align="center"><a href="../../README.md">简体中文</a> · <a href="README.en.md">English</a> · <strong>日本語</strong> · <a href="README.ko.md">한국어</a> · <a href="README.es.md">Español</a> · <a href="README.fr.md">Français</a> · <a href="README.de.md">Deutsch</a> · <a href="README.pt.md">Português</a> · <a href="README.ru.md">Русский</a> · <a href="README.ar.md">العربية</a></p>

OpenQuantum が完成させた6本の研究論文を収録しています。ランダム量子符号に共通する回復則から、Qwen の重みに組み込んだ量子回路シミュレータまで、結果、導出、検証記録を確認できます。

OpenQuantum を使って2量子ビット系の基底状態をローカルで再計算し、エネルギーを独立な厳密解と照合できます。モデルの API キーは不要です。

**[6本の論文を読む](#papers)** · **[2量子ビット系の基底状態を再計算][openquantum-first-run]** · [OpenQuantum をインストール][openquantum-start]

<a id="papers"></a>

## 6本の論文

論文と補足資料は英語の PDF です。

### QwenQPU · Qwen に量子回路シミュレータを組み込む

重みを固定した標準 Qwen レイヤー4層が、6量子ビット回路を古典的な状態ベクトル法でシミュレーションします。複素振幅を更新し、測定確率を読み出します。同じ言語モデルで実行機構を有効にすると、固定の48課題で基準を満たした数が **4/48 から 48/48** に増えました。Qiskit の参照値に対する確率の最大絶対誤差は **1.42 × 10⁻⁵** でした。

[論文と付録][qwenqpu] · [実装と評価データ][qwenqpu-code]

![QwenQPU のオリジナルバナー。回路命令、論理振幅の64個のアドレス、6量子ビットの古典シミュレーションで得た測定確率を示します。](../../assets/qwenqpu-banner.gif)

*右側の振幅と最終的な確率は、標準 Qwen レイヤーの実行記録に基づきます。*

### ランダム量子符号 · 古典情報と量子情報の回復に共通する臨界曲線

厳密なユニタリ 2-design から生成したランダム符号について、Pauli ノイズ下でメッセージ全体を復号する場合と、エンタングルメントを回復する場合の最適値を扱います。どちらの符号アンサンブル平均にも成り立つ有限符号長での上限と下限を導きました。所定の正則性を満たす独立同分布のノイズでは、正の符号化率における臨界領域で両者が**同じガウス曲線**に近づきます。

[論文][random-codes] · [補足資料][random-codes-si]

![青い帯が両方の課題の最適値のアンサンブル平均を挟みます。符号長が増すと帯は白い破線で示したガウス極限に向かって狭まります。](../../assets/random-code-banner.gif)

*論文の図 1b では、臨界領域の中心付近における帯の幅が、3種類の符号長にわたって約31.8から3.4パーセントポイントに縮まります。*

### 散逸する Majorana 系 · 四時刻相関でエネルギー緩和を取り出す

一定次数の相互作用と一様な散逸のもとで、4つのフェルミオン演算子の挿入を短時間の2組に分け、各組について全成分の和を取ります。2組の時間間隔がともに ε のとき、エネルギー信号は ε² のオーダーで、ほかの減衰モードの寄与は O(ε⁴) です。有限時間窓での誤差に上界を与えています。

[論文][majorana] · [補足資料][majorana-si]

![励起側と読み出し側でそれぞれフェルミオン対を全成分について足し合わせ、エネルギーモードを選び出します。右側は対の時間間隔を短くしたときの、ほかのモードの最大相対寄与です。](../../assets/majorana-banner.gif)

*図 1c の18個の有限系では、εJ = 0.16、gJ = 0.7 のとき、エネルギー信号に対するほかのモードの寄与は最大0.0983%でした。*

### フェルミオン系の Liouvillian 例外点 · 同じ初期状態から異なる減衰率

同じ初期状態でも、全粒子数とモード占有数の不均衡は異なる減衰チャネルを選びます。読み出しが保持された系または消失フラグに限られる場合、すべての状態準備を数えると、真空状態を用いた短時間試行が、ポンプ強度を識別するための最適な準備回数のスケーリングを達成します。

[論文][fermion-ep] · [補足資料][fermion-ep-si]

![同じ初期状態から、全粒子数とモード占有数の不均衡を別々に読み出すと、異なる減衰率が得られます。右側は条件付き信号と生存割合です。](../../assets/fermionic-ep-banner.gif)

*図 1 の12モードの例では、全粒子数とモード占有数の不均衡の緩和率はそれぞれ約1.679γ、1.118γです。γt = 5 での生存割合は 1.18 × 10⁻¹⁵ です。*

### リュードベリ鎖の E₇ 分光 · 2種類のプローブから4本の候補線

一様プローブと交互プローブを使い、偶パリティと奇パリティのチャネルを共通のエネルギー尺度で比較します。ブロッケードを厳密に課したモデルの特定のパラメータ範囲では、低エネルギー側の4つの候補ピークの質量比が E₇ の予測と **3% 以内**で一致します。36原子の計算点における浅い束縛状態 A₄ の結合エネルギーは理想値より約62%高く、有限鎖でのずれを定量化できます。

[論文][e7] · [補足資料][e7-si]

![交互プローブと一様プローブが、それぞれ2つの E₇ 候補ピークを捉えます。右側は36原子での両チャネルのスペクトルと浅い束縛状態の結合エネルギーです。](../../assets/rydberg-e7-banner.gif)

*図 2a と図 3a の36原子の計算では、A₄ 候補の結合エネルギーは最低エネルギーギャップの4.92%です。E₇ の理想値は3.04%です。*

### 開放リュードベリ系 · 減衰モードの重みの交換と密度信号

平均場の折り畳み点付近では、2つの減衰モードが回避交差を挟んで重みを交換します。コヒーレンスを消去することで、スペクトル方程式の根と観測量に対する重みを再現できます。準備条件を固定した48原子の例では、単一モードの重みの変化は約 **8.5 × 10⁻⁴**、密度信号全体の最大変化は約 **1.9 × 10⁻⁶** です。

[論文][rydberg-decay] · [補足資料][rydberg-decay-si]

![準備条件を固定すると、2つの減衰モードのスペクトル重みが交換されます。共通の対数目盛りで、同じ2つの離調点における単一モードの重みの変化と密度信号全体の最大差を示します。](../../assets/rydberg-decay-banner.gif)

*右側は図 5c の単一モードの重みの変化と図 5d の密度全体の最大差を、同じ密度の単位で比較しています。図中の2つの離調点の間で、変化量の大きさには約450倍の差があります。*

## OpenQuantum で自分の計算を実行する

OpenQuantum は量子研究向けのオープンソース AI プラットフォームです。Codex またはスタンドアロンのデスクトップアプリで使えます。量子回路を Qiskit でシミュレーションし、独立な方法で結果を検証して、同じ作業環境にあるコード、パラメータ、出力を次の計算に活用できます。

**[OpenQuantum をインストール][openquantum-start]** · [実行できる計算例を見る][openquantum-examples] · [GitHub の OpenQuantum][openquantum]

## 研究記録

[原稿の版と検証記録][manuscript-history] · [研究アーカイブ][research-archive] · [QwenQPU の原稿ソースと図][qwenqpu-tex]

[openquantum]: https://github.com/xi-zhao/OpenQuantum/blob/main/docs/readme/README.ja.md
[openquantum-start]: https://github.com/xi-zhao/OpenQuantum/blob/main/docs/readme/README.ja.md#quick-start
[openquantum-examples]: https://github.com/xi-zhao/OpenQuantum/blob/main/docs/readme/README.ja.md#algorithm-learning-and-method-comparison
[openquantum-first-run]: https://github.com/xi-zhao/OpenQuantum/blob/main/docs/readme/README.ja.md#no-model-key-needed-recompute-a-fixed-case
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
