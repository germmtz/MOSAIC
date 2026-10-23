# Predicting yield in a multi-neighbour grid

Use these grain-weight estimates for the same genotypes, at the model reference design conditions:

| Input file | Values to use |
| --- | --- |
| [DGE predictions](outputs/2026_IRD_experiment/ige_estimations/03_DGE_genotype_predictions.csv) | `prediction` = intercept + DGE; `prediction_SE` |
| [Pair interactions](outputs/2026_IRD_experiment/ige_estimations/03_all_genotype_interactions.csv) | `interaction` = neighbour IGE + focal-by-neighbour effect; `prediction_SE` |
| [Joint error covariance](outputs/2026_IRD_experiment/ige_estimations/prediction_error_covariance_poids_grains_g.csv) | Matrix **C**, matched by `prediction_id` in both rows and columns |

CSV files use semicolons and decimal commas (`read.csv2()` in R). Interactions are **directed**: `genotype` is the focal plant and `neighbour` is the affecting plant, including same-genotype pairs.

## Prediction

For focal genotype g, baseline B_g, neighbour genotypes h_j and chosen distance/rank weights w_j:

$$
\widehat Y_g=\widehat B_g+\sum_j w_j\widehat I_{g,h_j}.
$$

Illustration (weights shown on arrows):

```text
             117
              | 1
              v
  117 --1--> [35] <--0.5-- 233
```

Here, predicted yield is B_35 + 2 I_(35,117) + 0.5 I_(35,233). Add the baseline **once**; sum weights for repeated genotype pairs. The baseline comes from competition pots, not a solo-yield model. The current interactions contain no pot-type fixed contrast.

## Joint prediction uncertainty

**C is the joint prediction-error covariance matrix**, stored in
[`prediction_error_covariance_poids_grains_g.csv`](outputs/2026_IRD_experiment/ige_estimations/prediction_error_covariance_poids_grains_g.csv).
It describes uncertainty in the estimates from **both** `03_...` prediction files, including uncertainty shared between a baseline and an interaction or between two interactions. It is not a matrix of observed yield correlations.

Read the first column (`prediction_id`) as row labels, not numeric data. The remaining column headings are the matching prediction IDs. A diagonal entry is one estimate's `prediction_SE` squared (units: g²); an off-diagonal entry is the covariance between two estimation errors. Positive covariance increases uncertainty when the estimates are added; negative covariance reduces it.

**a is a vector of multipliers you build from the layout**, with one entry for every prediction ID in C. It contains no yields or SEs. For one focal plant, put 1 on its genotype baseline, the sum of neighbour weights on each applicable directed interaction, and 0 on every unused estimate. For the illustrated focal plant:

| Prediction ID (same labels in both CSVs and C) | Estimate in theta | Multiplier in a |
| --- | --- | ---: |
| `DGE:poids_grains_g:35` | `prediction` for genotype 35 | 1 |
| `interaction:poids_grains_g:35:117` | `interaction` for focal 35, neighbour 117 | 2 |
| `interaction:poids_grains_g:35:233` | `interaction` for focal 35, neighbour 233 | 0.5 |
| All other IDs | Their corresponding estimates | 0 |

Build theta by matching these IDs to `prediction` in the DGE file and `interaction` in the interactions file. **Order theta and a exactly like the rows and columns of C; do not assume the two prediction files have that order.** Then:

$$
\widehat Y=a^\top\widehat\theta,\qquad
SE(\widehat Y)=\sqrt{a^\top C a},\qquad
\text{approximate 95\% interval}=\widehat Y\pm1.96\,SE(\widehat Y).
$$

For a **whole layout**, construct one vector A by summing all focal plants' vectors a. Its baseline entries are the numbers of plants of each genotype; its interaction entries are the total directed neighbour weights for each focal–neighbour pair. Use A in the same formula for total yield and its SE. A neighbouring pair contributes in both directions when predicting both plants. For mean yield per plant, divide A by the number of plants; for a difference between layouts, use A_layout1 − A_layout2.

**Do not add individual SEs or confidence limits, or assume plants' prediction errors are independent.** Reusing an estimate across many plants reuses the same estimation error; C accounts for this shared uncertainty.

These intervals describe uncertainty in expected yield conditional on fitted variance parameters. They exclude future plant residual variation, model-selection uncertainty and extrapolation error. Neighbour weights, additivity and transfer to higher densities/third-genotype contexts require validation; use simulations as exploratory predictions.

