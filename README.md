# Dividing by a difference of noisy baselines

Code, data and software for *Dividing by a difference of noisy baselines: a
feasibility diagnostic for studies of behavioural coupling*.

Studies of behavioural coupling — mimicry, convergence, synchrony, accommodation
— observe two individuals together and compare against baselines estimated
separately for each. The coupling estimator therefore divides by a difference of
two noisy baselines, and this has three consequences that a conventional power
analysis does not surface:

- Precision is governed by a **solo contrast ratio** `nu`, not by the number of
  dyads. A thousand dyads of similar individuals carry no more information than
  ten, because the quantity being contracted is absent in each.
- If the observable is a count of events, its measurement error has an
  **irreducible floor** at `sqrt(total events)` that no instrument improves.
- If the grouping variable is an ordinal scale mapped proportionally to coupling
  strength, the attainable contrast in the observable is **capped** at
  `(sqrt(q) - 1) / (sqrt(q) + 1)` before any noise is considered.

These can combine to make a design unable to detect the effect *at all*, as
distinct from needing more data. This repository contains `couplingdx`, which
checks that before data collection, together with the analysis reproducing the
paper's worked case on within-meal bite mimicry.

## Contents

```
couplingdx/     the diagnostic, as an installable package (no dependencies)
analysis/       scripts reproducing the paper's fits, simulation and figures
figures/        Figures 1-4 as they appear in the paper
data/           derived bite timelines (see below)
```

## Install

```bash
cd couplingdx && pip install -e . && python -m pytest tests
```

No dependencies.

## Use

```python
from couplingdx import design_check

report = design_check(
    baseline_se=0.0158,        # SE of one person's baseline, one session
    baseline_contrast=0.0154,  # expected |difference| between the two
    gap=11.4,                  # expected joint-session gap under no coupling
    n_events_total=127,        # total events across both individuals
)
print(report)
```

```
Coupling design feasibility
==================================
  solo contrast ratio  nu   =     0.69   (need >= 12)
  expected gap              =    11.40
  gap measurement error     =    11.27  [counting floor]
  max spread across scale   =    0.451

  NOT FEASIBLE: the usable window in theta is empty.
```

Or from the command line:

```bash
python -m couplingdx --baseline-se 0.0158 --baseline-contrast 0.0154 \
                     --gap 11.4 --n-events-total 127
```

## How long a session would be needed?

For a counting observable the gap grows as `contrast * T` while its error grows
as `sqrt(rate * T)`, so the ratio grows as `sqrt(T)` and a sufficient duration
always exists in principle:

```python
from couplingdx import required_duration

required_duration(0.0154, 2*0.0908)                 # bite mimicry: 110,000 s (30 h)
required_duration(1.24, 2*5.12, dispersion=3.0)     # speech rate:    2,900 s (48 min)
```

Two designs of the same shape, differing by a factor of 38 in what they demand.
Contrast relative to the square root of the event rate is what decides
practicality — not either quantity alone.

## Choosing `nu_star`

`NU_MOMENTS = 3` is where the ratio estimator's moments begin to exist.
`NU_USABLE = 12` is where simulation shows recovery becomes acceptable
(>90% of replicates estimable, interquartile range below the true value).
They differ by a factor of four, and using the former as a design criterion is
the most common way to get this wrong. The default is `NU_USABLE`.

## Caveats

- `design_check` assumes independent baselines of equal precision. Correlated
  baselines (shared session, shared rater) inflate the variance further.
- An omitted `sigma` is taken as the counting floor, a lower bound, so the
  verdict is optimistic. Supply a measured value when you have one.
- `dispersion > 1` for clustered events (syllables in utterances, bites in
  mouthfuls); ignoring it understates the error.
- The thresholds come from simulation of one estimator. They are a guide, not a
  guarantee, and a simulation matched to your own design is better.

## Reproducing the paper

Requires `numpy`, `pandas`, `scipy`, `matplotlib`, `openpyxl`. Place the derived
timeline workbook at `data/modified bite_gt from FIC data (smm).xlsx`, then:

```bash
cd analysis
python poisson_mle_refit.py           # -> poisson_permeal.csv
python simulation_and_sensitivity.py  # simulation + sensitivity tables
python make_figures.py                # -> ../figures/fig1-4.pdf
```

`poisson_permeal.csv` and `permeal.csv` are committed, so the figures and
simulation can be reproduced without the workbook; only the first script needs it.

## Data

The derived timelines are built from the Food Intake Cycle dataset (Kyritsis,
Diou & Delopoulos, 2021). The underlying FIC recordings are the property of their
originators and are available from them under their own terms:
https://mug.ee.auth.gr/intake-cycle-detection/

Note that the derived timelines do **not** retain participant identifiers, which
is a known limitation discussed in the paper: no participant-level split,
mixed-effects fit, or exclusion of same-participant pseudo-dyads is possible from
these files alone.

## License

MIT for the code. The FIC data are under their originators' terms.
