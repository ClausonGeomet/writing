# Architecture — EMC feasibility writing and prototype

> Run `REQ-20261026-emc-part2`, 26 October 2026. This repository is a Quarto writing site; the research code is a bounded, standalone experiment, not a package API.

## System Architecture

### Module Structure

```mermaid
%%{init: {'theme': 'neutral'}}%%
graph TD
  subgraph Site["Writing site"]
    Q["_quarto.yml"]
    P["Published post"]
    D["Part 2 draft"]
  end
  subgraph Study["Study orchestration"]
    G["generation.R"]
    F["run_feasibility.py"]
    B["run_partb.py"]
    S["summarize.py"]
    V["provenance.py"]
  end
  subgraph Methods["Methods"]
    A["baselines.py"]
    O["emc_prototypes.py"]
  end
  subgraph Evidence["Frozen evidence"]
    X["fixtures and C"]
    Y["A results and freeze"]
    Z["B results and manifest"]
  end
  G --> X
  F --> X
  F --> A
  F --> Y
  B --> X
  B --> Y
  B --> O
  B --> Z
  S --> Y
  V --> Y
  V --> Z
  D --> Y
  D --> Z
  Q --> P
  style D fill:#1e90ff,stroke:#1565c0,color:#fff
```

| Module / file | Purpose | Key dependencies | Changed in this run |
| --- | --- | --- | --- |
| `_quarto.yml`, published post | Site render list and Part 1 context | Quarto | No |
| `drafts/emc-part2-review.qmd` | Review article excluded from site render | Frozen A/B evidence | Yes, docs phase |
| `research/emc-part2a/generation.R` | Known-truth ideal assay generation | Pinned simCoDaR | Yes, simulator phase |
| `research/emc-part2a/run_feasibility.py`, `baselines.py` | Freeze fixtures and five baseline fits | NumPy, SciPy | Yes, simulator phase |
| `research/emc-part2a/run_partb.py` | B fit orchestration and held-out scoring | Frozen A, proposed solver | Yes, simulator phase |
| `research/emc-part2a/summarize.py`, `reporting_provenance.py` | Aggregate metrics and provenance derivatives | CSVs, freeze manifest | Yes, simulator phase |
| `research/emc-part2b/prototypes/emc_prototypes.py` | Allocation, scoring, finite order search | NumPy | Yes, builder phase |
| `research/emc-part2a/fixtures`, `results`; `research/emc-part2b/results` | Input, predictions, hashes and findings | Prior stages | Yes, simulator/tester phase |

### Function Call Graph

```mermaid
%%{init: {'theme': 'neutral'}}%%
graph TD
  Gen["simulate scenario"] --> Conv["minerals_to_assay"]
  RunA["run feasibility"] --> Solve["solve baseline"]
  RunA --> Rank["rank and LP bounds"]
  RunB["run Part B"] --> Fit["fit_order"]
  Fit --> Alloc["allocate"]
  Fit --> Score["score_assays"]
  RunB --> Predict["predict_order"]
  Predict --> Alloc
  RunB --> Metrics["held-out metrics"]
  style Fit fill:#1e90ff,stroke:#1565c0,color:#fff
  style Alloc fill:#1e90ff,stroke:#1565c0,color:#fff
  style Score fill:#1e90ff,stroke:#1565c0,color:#fff
  style Predict fill:#1e90ff,stroke:#1565c0,color:#fff
```

| Function | Defined in | Called by | Purpose | Changed |
| --- | --- | --- | --- | --- |
| `simulate_mineral_assay_scenario`, `minerals_to_assay` | simCoDaR dependency | `generation.R` | Generate truth and oxide response | No |
| `solve` | `baselines.py` | `run_feasibility.py` | Five rowwise baseline fits | Yes |
| `fit_order` | `emc_prototypes.py` | `run_partb.py` | Train global phase order from assays | Yes |
| `allocate`, `predict_order` | `emc_prototypes.py` | `fit_order`, `run_partb.py` | Limiting component fit and replay | Yes |
| `score_assays` | `emc_prototypes.py` | `fit_order`, `run_partb.py` | Row and distribution losses | Yes |

### Data Flow

```mermaid
%%{init: {'theme': 'neutral'}}%%
graph TD
  Source["Pinned simCoDaR"] --> Fixture["Truth and ideal assays"]
  Fixture --> Split["30 / 10 / 20 split"]
  Split --> Baseline["Fit and freeze A"]
  Baseline --> Check{{"Freeze verified?"}}
  Check -- yes --> Search["Train B order"]
  Check -- no --> Stop["Stop evaluation"]
  Search --> Heldout["Score held-out rows"]
  Heldout --> Report["Provenance and prose"]
  style Report fill:#1e90ff,stroke:#1565c0,color:#fff
```

| Stage | Input | Output | Changed |
| --- | --- | --- | --- |
| Fixture generation | Pinned R source and seeds | Truth, assay rows, response matrices | Yes |
| A freeze | Fixtures and five solvers | Frozen baseline results and hash ledger | Yes |
| B evaluation | Verified A freeze and training assays | Orders, predictions, metrics, traces | Yes |
| Documentation | Audited CSVs and source ledger | Research reports and excluded Quarto draft | Yes |

## Notes

Part B fitting sees assays and chemistry, not mineral truth; truth enters held-out evaluation. The draft remains outside the site's render list. The historical v2 snapshot has one unavailable referenced manifest (55/56); the current freeze verifies 55/55.
