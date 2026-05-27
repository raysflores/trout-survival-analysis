# Survival Analysis of Stocked Brook and Brown Trout in Pennsylvania Streams

**Raymond Santiago Flores · STAT 3559: Survival Analysis · Spring 2026 · University of Virginia**

A time-to-event analysis of post-stocking survival in hatchery-reared trout across five central Pennsylvania streams, applying Kaplan-Meier estimation, log-rank testing, and Cox Proportional Hazards regression to identify behavioral and environmental predictors of mortality risk.

---

## Overview

State fish stocking programs release millions of hatchery-reared trout annually, yet post-stocking survival rates remain poorly understood. Fish that disperse rapidly from stocking sites face elevated predation and displacement risk; fish that hold position may benefit from localized resources but remain vulnerable to angling. This analysis uses telemetry-derived survival data from a U.S. Fish and Wildlife Service / data.gov dataset to model time-to-mortality for 90 individually tracked trout across five streams in Pennsylvania's Ridge and Valley physiographic province.

The study period spans 153 days (maximum observed survival: 166 days). Of 90 fish in the analytic sample (99 original observations minus 9 excluded due to missing `move_rate` data), **62 experienced the mortality event (68.9% event rate)** and 28 were right-censored at the study end date or upon confirmed emigration beyond the monitored reach.

---

## Data

**Source:** U.S. Fish and Wildlife Service — publicly available via data.gov

**Analytic sample:** n = 90 (after exclusion of 9 observations missing `move_rate`)

| Variable | Description |
|---|---|
| Survival time (days) | Days from stocking to mortality or censoring |
| Event indicator | 1 = mortality observed; 0 = right-censored |
| `species` | Brook trout (*Salvelinus fontinalis*, n = 47) or Brown trout (*Salmo trutta*, n = 43) |
| `Wgt` | Individual fish weight (g) at stocking; log-transformed in models |
| `Lgt` | Fork length (mm) at stocking — excluded due to high collinearity with `Wgt` (r = 0.935) |
| `move_rate` | Mean daily movement rate (m/day) — log-transformed in models |
| `order` | Stream order category: low (n = 7), medium (n = 54), complex (n = 29) |
| `cumdrain` | Cumulative drainage area of the stocking reach (km²) |
| `stream` | Individual stream identity: Hunts Run (n = 29), McKinnon Run (n = 18), McNuff Branch (n = 18), Rock Run (n = 7), Whitehead Run (n = 18) |

> **Note:** `stream` and `order` are perfectly collinear — each stream maps to exactly one order category — and were not included in the same model. `order` was retained as the theoretically preferred grouping variable.

### Survival Time Summary

| Statistic | Value |
|---|---|
| Minimum | 16 days |
| Q1 | 43 days |
| Median | 98 days |
| Mean | 98.2 days |
| Q3 | 153 days |
| Maximum | 166 days |
| SD | 51.3 days |

---

## Methods

### Kaplan-Meier Estimation and Log-Rank Test

Kaplan-Meier survival curves were estimated separately for Brook and Brown trout. A log-rank test assessed whether the species-specific survival functions differed significantly over the study period.

**Result:** χ² = 0.5, df = 1, **p = 0.5** — no statistically significant difference in survival between species.

### Cox Proportional Hazards Regression

Three nested Cox PH models were fit and compared by AIC. Continuous predictors `Wgt` and `move_rate` were log-transformed to satisfy linearity-in-log-hazard assumptions. Fork length (`Lgt`) was excluded from all models due to its high correlation with `Wgt` (r = 0.935), which would induce severe multicollinearity.

| Model | Variables | AIC |
|---|---|---|
| M1 | `species + log_Wgt + order + cumdrain` | 512.80 |
| M2 | `species + log_Wgt + stream + cumdrain` | 512.80 |
| **M3** (final) | `species + log_Wgt + order + cumdrain + log_move_rate` | **505.87** |

Model 3 was selected as the final specification based on lowest AIC (Δ AIC = 6.93 vs. M1/M2) and concordance index C = 0.60.

### Final Model (M3) Results

```
Surv(time, event) ~ species + log_Wgt + order + cumdrain + log_move_rate
```

| Predictor | Hazard Ratio | 95% CI | p-value |
|---|---|---|---|
| species (Brown vs. Brook) | — | — | ns |
| log_Wgt | — | — | ns |
| order (medium vs. low) | — | — | ns |
| order (complex vs. low) | — | — | ns |
| cumdrain | — | — | ns |
| **log_move_rate** | **0.793** | — | **0.032** |

`log_move_rate` was the only statistically significant predictor (p = 0.032). A one-unit increase in log(move_rate) is associated with a **21% decrease in the instantaneous hazard of mortality** (HR = 0.793), suggesting that fish with higher daily movement rates have improved survival — potentially reflecting active habitat selection behavior or faster emigration from high-risk areas.

### Proportional Hazards Assumption

The proportional hazards assumption was evaluated using Schoenfeld residuals for each covariate and globally. The PH assumption was **violated for all key predictors:**

| Covariate | Test p-value |
|---|---|
| `species` | 0.0009 |
| `order` | 0.016 |
| `cumdrain` | 0.044 |
| `log_move_rate` | 0.0015 |
| **Global test** | **6.5 × 10⁻⁵** |

These violations indicate that hazard ratios change over the study period and are not constant — a fundamental violation of the standard Cox model. The study journal discusses two methodological corrections appropriate for future analysis:

1. **Time-varying covariates:** Allow covariate effects to vary as a function of time by including time × covariate interaction terms
2. **Stratified Cox model:** Stratify on the variable(s) most severely violating PH (e.g., `species`, `order`) to allow stratum-specific baseline hazard functions while retaining shared coefficient estimates for other predictors

The PH violations represent the primary limitation of the current analysis and are a focus of the interpretive discussion.

---

## Technical Stack

- **Language:** R
- **Key packages:** `survival`, `survminer`, `dplyr`, `tidyr`, `ggplot2`, `knitr`, `kableExtra`, `gridExtra`
- **Methods:** Kaplan-Meier estimation, log-rank test, Cox Proportional Hazards regression, AIC model comparison, Schoenfeld residual test (PH assumption), concordance index, right-censoring

---

## Repository Structure

```
├── Project Script.Rmd         # Full analysis — KM curves, Cox models, PH diagnostics
├── data/
│   └── troutsurvival_daily.csv   # Individual-level telemetry data
├── deliverables/
│   └── TroutFinalJournal.pdf     # Final written analysis journal (22 pp.)
└── README.md
```

---

## Study Site

Five streams in the **Ridge and Valley physiographic province, central Pennsylvania**:

| Stream | n | Stream Order |
|---|---|---|
| Hunts Run | 29 | — |
| McKinnon Run | 18 | — |
| McNuff Branch | 18 | — |
| Rock Run | 7 | — |
| Whitehead Run | 18 | — |

Stream order categories (low / medium / complex) are defined in the source dataset and map one-to-one to individual streams.

---

## Limitations

- **PH assumption violations:** Hazard ratios are non-constant over time for all modeled predictors (global Schoenfeld p = 6.5 × 10⁻⁵). Standard Cox estimates should be interpreted as average effects rather than constant hazard ratios. Time-varying or stratified extensions are recommended.
- **Small n in some strata:** Rock Run (n = 7) has insufficient observations for stable stratum-specific estimation.
- **Censoring mechanism:** Right-censoring assumes non-informative censoring (i.e., censored fish do not differ systematically from uncensored fish in their underlying survival distribution). Emigration-based censoring may violate this assumption if movement-prone fish are both more likely to emigrate and to survive.
- **Single stocking cohort:** Results may not generalize across years, stocking densities, or stream conditions outside the study period.

---

## About

Completed as the course project for STAT 3559: Survival Analysis at the University of Virginia, Spring 2026. Data sourced from the U.S. Fish and Wildlife Service via data.gov.

**Author:** Raymond Santiago Flores | [raymondosf@gmail.com](mailto:raymondosf@gmail.com)
