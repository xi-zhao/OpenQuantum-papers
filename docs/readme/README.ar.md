# OpenQuantum Papers

<p align="center"><a href="../../README.md">简体中文</a> · <a href="README.en.md">English</a> · <a href="README.ja.md">日本語</a> · <a href="README.ko.md">한국어</a> · <a href="README.es.md">Español</a> · <a href="README.fr.md">Français</a> · <a href="README.de.md">Deutsch</a> · <a href="README.pt.md">Português</a> · <a href="README.ru.md">Русский</a> · <strong>العربية</strong></p>

ست أوراق بحثية أنجزها OpenQuantum، من قانون استعادة مشترك للشفرات الكمومية العشوائية إلى محاكٍ للدوائر الكمومية مدمج في أوزان Qwen. تجد هنا النتائج والاشتقاقات وسجلات التحقق.

أعد حساب الحالة القاعية لنظام من كيوبتين محليًا باستخدام OpenQuantum، وقارن الطاقة بحل دقيق مستقل. يعمل المثال من دون مفتاح للنموذج.

**[اقرأ الأوراق الست](#papers)** · **[أعد حساب الحالة القاعية لكيوبتين][openquantum-first-run]** · [ثبّت OpenQuantum][openquantum-start]

<a id="papers"></a>

## الأوراق الست

تتوفر الأوراق وملاحقها باللغة الإنجليزية بصيغة PDF.

### محاكاة الدوائر الكمومية داخل Qwen · QwenQPU

تنفذ أربع طبقات قياسية من Qwen بأوزان ثابتة محاكاة كلاسيكية بمتجه الحالة لدوائر من ستة كيوبتات. تحدّث الطبقات السعات المركبة وتقرأ احتمالات القياس. عند تفعيل آلية التنفيذ في نموذج اللغة نفسه، ارتفع عدد المهام التي حققت معيار النجاح في التقييم الثابت من <strong dir="ltr">4/48</strong> إلى <strong dir="ltr">48/48</strong>. بلغ أكبر خطأ مطلق في الاحتمالات مقارنةً بمرجع Qiskit نحو <strong dir="ltr">1.42 × 10⁻⁵</strong>.

[الورقة والملاحق][qwenqpu] · [التنفيذ وبيانات التقييم][qwenqpu-code]

![الصورة المتحركة الأصلية لـ QwenQPU. تعرض تعليمات الدائرة و64 عنوانًا منطقيًا للسعات واحتمالات ناتجة من محاكاة كلاسيكية لستة كيوبتات.](../../assets/qwenqpu-banner.gif)

*السعات والاحتمال النهائي في الجزء الأيمن مأخوذة من تشغيل فعلي لطبقات Qwen القياسية.*

### الشفرات الكمومية العشوائية · منحنى حرج مشترك لاستعادة المعلومات الكلاسيكية والكمومية

تشتق الورقة حدودًا علوية وسفلية مشتركة عند طول كتلة محدود لشفرات عشوائية مسحوبة من تصميم وحدوي دقيق من الرتبة الثانية. تنطبق هذه الحدود على متوسط الأداء الأمثل عبر مجموعة الشفرات في كلٍّ من فك ترميز الرسالة كاملةً واستعادة التشابك تحت ضجيج باولي. وعندما يكون الضجيج مستقلًا ومتماثل التوزيع ويحقق شروط الانتظام، يقترب المتوسطان من **المنحنى الغاوسي نفسه** في النافذة الحرجة عند معدل ترميز موجب.

[الورقة][random-codes] · [المواد التكميلية][random-codes-si]

![يحصر الشريط الأزرق متوسطَي الأداء الأمثل للمهمتين. ومع ازدياد طول الكتلة، يضيق الشريط باتجاه الحد الغاوسي المرسوم بخط أبيض متقطع.](../../assets/random-code-banner.gif)

*قرب مركز النافذة الحرجة، ينخفض عرض الشريط من نحو 31.8 إلى 3.4 نقطة مئوية عبر أطوال الكتل الثلاثة في الشكل 1b.*

### أنظمة ماجورانا التبددية · ارتباط رباعي الأزمنة ينتقي استرخاء الطاقة

في ظل تآثرات من رتبة ثابتة وتبدد متجانس، تُجمع الإدخالات الفرميونية الأربعة في زوجين بفاصلين زمنيين قصيرين، ويُجمع كل زوج على جميع المركبات. عندما يساوي الفاصل في كل زوج ε، تكون إشارة الطاقة من رتبة <span dir="ltr">ε²</span>، وتكون مساهمة أنماط الاضمحلال الأخرى من رتبة <span dir="ltr">O(ε⁴)</span>. تعطي الورقة حدًا للخطأ في نافذة زمنية محدودة.

[الورقة][majorana] · [المواد التكميلية][majorana-si]

![يؤدي جمع الزوج الفرميوني عند الإثارة والزوج عند القراءة على جميع المركبات إلى انتقاء نمط الطاقة. يعرض الجزء الأيمن أكبر مساهمة نسبية للأنماط الأخرى عند تقصير الفاصل الزمني.](../../assets/majorana-banner.gif)

*في الأنظمة المحدودة الـ18 في الشكل 1c، عند <span dir="ltr">εJ = 0.16</span> و<span dir="ltr">gJ = 0.7</span>، بلغت أقصى مساهمة للأنماط الأخرى نسبةً إلى إشارة الطاقة <span dir="ltr">0.0983%</span>.*

### النقاط الاستثنائية لمؤثر ليوفيلي فرميوني · معدّلا اضمحلال من حالة ابتدائية واحدة

ينتقي العدد الكلي للجسيمات واختلال إشغال الأنماط قناتين مختلفتين للاضمحلال انطلاقًا من الحالة الابتدائية نفسها. إذا كانت القراءة تعيد النظام المحتفَظ به أو إشارة محو، فإن التجارب القصيرة التي تبدأ من حالة الفراغ تحقق التحجيم الأمثل لعدد التحضيرات اللازمة للتمييز بين شدّات الضخ، مع احتساب كل عملية تحضير.

[الورقة][fermion-ep] · [المواد التكميلية][fermion-ep-si]

![يعطي العدد الكلي للجسيمات واختلال إشغال الأنماط معدّلي اضمحلال مختلفين من الحالة الابتدائية نفسها. ويعرض الجزء الأيمن الإشارة المشروطة ونسبة البقاء.](../../assets/fermionic-ep-banner.gif)

*في مثال الأنماط الـ12 في الشكل 1، يبلغ معدّلا استرخاء العدد الكلي للجسيمات واختلال الإشغال نحو <span dir="ltr">1.679γ</span> و<span dir="ltr">1.118γ</span> على الترتيب. وعند <span dir="ltr">γt = 5</span> تبلغ نسبة البقاء <span dir="ltr">1.18 × 10⁻¹⁵</span>.*

### مطيافية E₇ في سلاسل ريدبرغ · مجسّان يكشفان أربعة خطوط مرشحة

يضع المجسّان المنتظم والمتناوب قناتي التكافؤ الزوجي والفردي على مقياس طاقة واحد. ضمن نطاق من معلمات نموذج الحجب الصارم، تختلف نسب الكتل لأربع قمم مرشحة منخفضة الطاقة عن تنبؤات E₇ بما لا يتجاوز <strong dir="ltr">3%</strong>. عند نقطة التشغيل ذات 36 ذرة، تزيد طاقة ربط المرشح ضعيف الارتباط A₄ بنحو <span dir="ltr">62%</span> على القيمة المثالية، وهو مقياس كمي لانحراف السلسلة المحدودة.

[الورقة][e7] · [المواد التكميلية][e7-si]

![يكشف كلٌّ من المجسّين المتناوب والمنتظم عن قمتين مرشحتين من طيف E₇. ويعرض الجزء الأيمن الطيفين وطاقة ربط الحالة ضعيفة الارتباط في نظام من 36 ذرة.](../../assets/rydberg-e7-banner.gif)

*في حساب الذرات الـ36 في الشكلين 2a و3a، تساوي طاقة ربط المرشح A₄ نسبة <span dir="ltr">4.92%</span> من أصغر فجوة طاقة. وتبلغ قيمة E₇ المثالية <span dir="ltr">3.04%</span>.*

### تجمّعات ريدبرغ المفتوحة · من تبادل أوزان أنماط الاضمحلال إلى إشارة الكثافة

قرب نقطة الطي في حل المجال المتوسط، يتبادل نمطان للاضمحلال وزنيهما عند تقاطع متجنّب. وباستبعاد مركبات التماسك حسابيًا، يمكن استعادة الجذور الطيفية والأوزان في الكمية المرصودة. في المثال المؤلف من 48 ذرة مع تثبيت التحضير، يتغير وزن نمط واحد بنحو <strong dir="ltr">8.5 × 10⁻⁴</strong>، ويبلغ أكبر تغير في إشارة الكثافة الكاملة نحو <strong dir="ltr">1.9 × 10⁻⁶</strong>.

[الورقة][rydberg-decay] · [المواد التكميلية][rydberg-decay-si]

![يتبادل نمطان للاضمحلال الوزن الطيفي مع تثبيت التحضير. يقارن مقياس لوغاريتمي مشترك تغير وزن أحد النمطين بأكبر فرق في إشارة الكثافة الكاملة بين القيمتين نفسيهما للانحراف عن الرنين.](../../assets/rydberg-decay-banner.gif)

*يقارن الجزء الأيمن تغير وزن النمط في الشكل 5c بأكبر فرق في الكثافة في الشكل 5d، بوحدة الكثافة نفسها. وبين قيمتي الانحراف المحددتين يبلغ التغير الأول نحو 450 ضعف الثاني.*

## نفّذ حساباتك باستخدام OpenQuantum

OpenQuantum منصة ذكاء اصطناعي مفتوحة المصدر للبحث الكمومي، تتوفر في Codex وكتطبيق مستقل لسطح المكتب. حاكِ دائرة كمومية باستخدام Qiskit، وتحقق من النتيجة بطريقة مستقلة، ثم واصل العمل بالشيفرة والمعلمات والمخرجات نفسها في مساحة عملك.

**[ثبّت OpenQuantum][openquantum-start]** · [تصفّح الأمثلة القابلة للتشغيل][openquantum-examples] · [OpenQuantum على GitHub][openquantum]

## سجلات البحث

[إصدارات المخطوطات وسجلات التحقق][manuscript-history] · [أرشيف البحث][research-archive] · [مصادر مخطوطة QwenQPU وأشكالها][qwenqpu-tex]

[openquantum]: https://github.com/xi-zhao/OpenQuantum/blob/main/docs/readme/README.ar.md
[openquantum-start]: https://github.com/xi-zhao/OpenQuantum/blob/main/docs/readme/README.ar.md#quick-start
[openquantum-examples]: https://github.com/xi-zhao/OpenQuantum/blob/main/docs/readme/README.ar.md#algorithm-learning-and-method-comparison
[openquantum-first-run]: https://github.com/xi-zhao/OpenQuantum/blob/main/docs/readme/README.ar.md#no-model-key-needed-recompute-a-fixed-case
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
