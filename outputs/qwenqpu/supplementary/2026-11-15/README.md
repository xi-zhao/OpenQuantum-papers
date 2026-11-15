# QwenQPU Supplementary Data

Supplement to **QwenQPU: A Large Language Model with an Internal Quantum Processing Unit**, by Xi Zhao and Yunlong Yu.

- Xi Zhao: CAS Key Laboratory of Quantum Information, University of Science and Technology of China
- Yunlong Yu: State Key Laboratory of Low Dimensional Quantum Physics, Department of Physics, Tsinghua University

## Download and verify

[Download the complete Supplementary Data](QwenQPU-engineering-supplement-20261115.zip?raw=1) (14,900,920 bytes). The package SHA-256 is:

```text
73b978119bf382cd1cc33394f8a7960abd49a1486f906acae1de808ebfa5c1a2
```

The archive contains the matched-executor and full-state records, construction tests and per-gate states, profiling data, exact experiment scripts, and an offline reconstruction entry point. It also includes the paired normalization analysis and an execution-provenance note. Model weights are supplied separately through the pinned sources described in the reproduction guide.

## Reproduce the saved results

Extract the ZIP, enter its `QwenQPU-engineering-supplement-20261115` directory, and run with Python 3.11 or later:

```bash
python3 reproduce.py --out ./reproduced-check
```

Use a new output directory. The verifier uses the Python standard library and requires no network, GPU, credentials, or model inference. It checks the evidence hashes and independently reconstructs the supplied state and readout results:

| Material | Coverage |
|---|---:|
| Matched requests and complete-model states | 48 tasks |
| Standalone QPU construction tests | 36 cases |
| Initialization and per-gate state snapshots | 927 snapshots |
| Paired uncorrected/corrected readouts | 96 events from 48 tasks |

The extracted `README.md` distinguishes offline verification, the standard-layer CPU check, complete-model demonstrations, and the archived GPU study with its original environment and directory dependencies. `PROVENANCE.md` documents the harness correction and reused control records. These records do not add independent samples to the benchmark.

## Data identity

The nested evidence archive is preserved byte for byte with SHA-256 `f3b9d8ea845cb160ac9470f44fb0e63b5e9b338b40bc30f33c6a8b497f09d963`. Its historical filenames, attribution, and operational notes are retained as provenance. The current article's authors are Xi Zhao and Yunlong Yu.

The execution-access ablation is separately identified by the r5 records at [QwenQPU commit caae45327b46458163d0840a4c1466e5575a7b4f](https://github.com/xi-zhao/QwenQPU/tree/caae45327b46458163d0840a4c1466e5575a7b4f). Its timing scope differs from the matched internal/external experiment supplied here.

`MANIFEST.json` records the package identity and authors; `SHA256SUMS` identifies this directory's files. Use the commit-specific URL in the article when citing these data.
