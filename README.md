# Predicting death in primary biliary cholangitis: prognostic factors and the effect of D-penicillamine

## Abstract

Primary biliary cholangitis (PBC) is a chronic autoimmune liver disease that slowly progresses to
cirrhosis. We used data from the Mayo Clinic PBC study (1974 to 1984; 418 patients, 125 deaths in
the randomised subset) to ask which baseline measurements predict death, and whether D-penicillamine
improves survival. Kaplan-Meier curves showed strong differences by histologic stage (log-rank
p < 0.001) but none by treatment (p = 0.75). In a multivariable Cox model, higher bilirubin (hazard
ratio 2.26 per log unit), higher protime and oedema were associated with death, while albumin was
protective. D-penicillamine had no clear effect (HR 0.89, 95% CI 0.62 to 1.29). A model using all
candidate variables predicted death in held-out patients with concordance 0.86, compared with 0.82
for bilirubin alone.

## 1. Introduction

PBC destroys the small bile ducts of the liver. It mostly affects middle-aged women, and it progresses
slowly until cirrhosis and liver failure develop. Clinicians need to know which patients are at high risk
of death, so that they can plan monitoring and decide when to refer for transplantation. Before the
introduction of ursodeoxycholic acid, D-penicillamine was tested as a treatment.

Our aim was to answer a data-driven question: **which baseline clinical and laboratory measurements
predict death in PBC, and does D-penicillamine change survival?** We expected serum bilirubin and albumin
to be the strongest markers of liver function, and we expected histologic stage to matter. We tested
whether treatment has any effect in the randomised trial.

## 2. Data

**Source.** The `pbc` dataset in the `survival` R package (version 3.5-8; Therneau, 2024). It comes from
the Mayo Clinic trial of D-penicillamine against placebo, described in Fleming and Harrington (1991) and
Therneau and Grambsch (2000). The data are exported to `data/pbc.csv` by `code/00_export_data.R`.

**Size and variables.** 418 patients, registered 1974 to 1984. The first 312 were randomised to
D-penicillamine (158) or placebo (154). The other 106 were not randomised and have fewer measurements.
Outcome variables are `time` (days from registration to death, transplant or the July 1986 analysis date)
and `status` (0 censored, 1 transplant, 2 death). There were 161 deaths, 25 transplants and 232 patients
censored alive across the full dataset. Clinical variables include age, sex, ascites, hepatomegaly,
spider angiomata, oedema, serum bilirubin, albumin, prothrombin time (protime) and histologic stage (1 to 4).

**Quality.** Missing values are concentrated in the 106 non-randomised patients, who lack treatment, ascites,
hepatomegaly, spiders, alk.phos and AST. Triglycerides, cholesterol and copper are missing for 26% to 33%
of patients (`results/missingness.csv`). The 312 randomised
patients have complete data for every variable in the main model. Stage is missing for 6 patients, and these
are excluded from the stage analysis only.

**Ethics and licence.** The data are de-identified and included in a public R package. The `survival`
package is licensed under LGPL (version 2 or later). The original trial consent is not described in the
package documentation, so we rely on the package maintainers' release. No new patient data were collected.

## 3. Methods

The workflow is in `code/00_export_data.R` (export) and `code/01_analysis.R` (analysis).

**Import and QC.** We read the CSV, counted outcomes and missing values, and checked that the coding
matched the package documentation: treatment 1 is D-penicillamine and 2 is placebo, and status 2 is death.

**Outcome definition.** The event is death. Transplantation is treated as censoring at the time of transplant,
since transplanted patients are no longer followed for death. This is the usual approach for this dataset,
but it is a choice, and it is discussed below.

**Exploratory analysis.** We plotted follow-up time by outcome (Fig. 1) and bilirubin by outcome on a log
scale (Fig. 2). Bilirubin is strongly right-skewed, so we used its logarithm in the models. Baseline
characteristics were compared by treatment arm (`results/baseline_by_treatment.csv`).

**Survival curves and tests.** Kaplan-Meier estimates were drawn overall, by histologic stage and by randomised
treatment (Figs 3 and 4). The log-rank test compared groups.

**Cox proportional hazards models.** Each variable was first tested alone (univariable). A multivariable model
then included treatment, age, sex, ascites, hepatomegaly, spiders, oedema, log bilirubin, albumin, log protime
and stage. It was fitted to the 312 randomised patients with complete data. The proportional hazards assumption
was checked with Schoenfeld residuals (`cox.zph`; Fig. 6).

**Validation.** Patients were split 70/30, stratified by death, to give 219 training and 93 test patients
(37 deaths in the test set). The multivariable model and a bilirubin-only model were fitted on the training set.
Both were scored on the test set with Harrell's concordance index (C; Harrell et al., 1996). A C of 0.5 means
chance and 1 means perfect ranking of risk.

**R packages.** `survival` (Therneau, 2024) for survival analysis, `tidyverse` (Wickham et al., 2019) for data
handling, and `ggplot2` (Wickham, 2016) for graphics. Versions are in `results/session_info.txt`.

## 4. Results

### 4.1 Patients and follow-up

Of the 418 patients, 125 died among the 312 randomised patients, and 161 died in the full dataset. Median
follow-up was 4.7 years in the whole dataset (`results/qc_summary.csv`). Overall median survival, estimated by
Kaplan-Meier, was 9.3 years.



<img width="1050" height="675" alt="fig1_followup" src="https://github.com/user-attachments/assets/510870b1-937b-4b6a-b2c3-a3828ebc4e8f" />


*Fig 1. Follow-up time by outcome, across all 418 patients.*


Bilirubin was higher in patients who died (Fig. 2).


<img width="900" height="675" alt="fig2_bilirubin_by_outcome" src="https://github.com/user-attachments/assets/7e3c6260-a09d-4c03-ae9f-47251e6b99d4" />


*Fig 2. Serum bilirubin by outcome, on a log scale.*

Treatment groups were similar at baseline (Table 1).

| Arm | Patients | Deaths | Mean age (years) | Female (%) | Median bilirubin (mg/dl) | Median albumin (g/dl) |
|---|---|---|---|---|---|---|
| Placebo | 154 | 60 | 48.6 | 90.3 | 1.3 | 3.54 |
| D-penicillamine | 158 | 65 | 51.4 | 86.7 | 1.4 | 3.56 |

*Table 1. Baseline characteristics of the randomised subset (`results/baseline_by_treatment.csv`).*

### 4.2 Survival by stage and treatment

Survival differed strongly by histologic stage (log-rank p = 4.1 x 10^-15; Fig. 3). Stage 4 patients had the
lowest survival, and stage 1 patients the highest. Stage 1 has only 21 patients, so its curve is imprecise.


<img width="1050" height="750" alt="fig3_km_stage" src="https://github.com/user-attachments/assets/763da72b-cd59-4f23-bbe2-3824717313c4" />


*Fig 3. Kaplan-Meier survival by histologic stage, with 95% confidence bands. Log-rank p = 4.1 x 10^-15.*

Survival did not differ by treatment (log-rank p = 0.75; Fig. 4).


<img width="1050" height="750" alt="fig4_km_treatment" src="https://github.com/user-attachments/assets/443940ff-25d4-4fbb-a70b-f1c192a4abdf" />


*Fig 4. Kaplan-Meier survival by randomised treatment (312 patients). Log-rank p = 0.75.*

### 4.3 Factors associated with death

In the multivariable model (Fig. 5, `results/cox_multivariable.csv`), the following were associated with death:

- **Bilirubin**: HR 2.26 per unit of log bilirubin (95% CI 1.81 to 2.81; p < 0.001).
- **Albumin**: HR 0.47 per g/dl (95% CI 0.27 to 0.79; p = 0.005). Higher albumin was protective.
- **Age**: HR 1.03 per year (95% CI 1.01 to 1.05; p = 0.007).
- **Oedema**: HR 2.19 per level (95% CI 1.19 to 4.03; p = 0.012).
- **Protime**: HR 20.0 per unit of log protime (95% CI 1.84 to 217.7; p = 0.014). This estimate is very imprecise.

Treatment with D-penicillamine had HR 0.89 (95% CI 0.62 to 1.29; p = 0.54). Stage was not significant
overall, though stage 4 showed a trend. Sex, ascites, hepatomegaly and spiders were not significant.


<img width="1050" height="825" alt="fig5_forest" src="https://github.com/user-attachments/assets/6c1fa64a-18c1-418a-822f-477289d854b3" />


*Fig 5. Adjusted hazard ratios from the multivariable Cox model, on a log scale. Stage coefficients are
relative to stage 1.*

**Proportional hazards.** The global test was not significant (p = 0.15). Two individual terms showed departures
from proportional hazards: oedema (p = 0.028) and log protime (p = 0.018). The hazard ratios for these two
variables are averages over follow-up and should be read with care. Fig. 6 shows the Schoenfeld residuals for
bilirubin, which met the assumption (p = 0.34).


<img width="1000" height="700" alt="fig6_schoenfeld_log_bili" src="https://github.com/user-attachments/assets/63130038-7d70-46d7-97b8-ca4e803658d5" />


*Fig 6. Schoenfeld residuals for log bilirubin. A flat trend supports the proportional hazards assumption.*

### 4.4 Prediction

On held-out patients (Table 2), the multivariable model had concordance 0.86. Bilirubin alone had 0.82.

| Model | C (training) | C (held-out test) |
|---|---|---|
| Multivariable Cox (all candidate variables) | 0.846 | 0.863 |
| Bilirubin only | 0.785 | 0.821 |

*Table 2. Concordance on 219 training and 93 held-out patients (37 deaths in the test set).
From `results/model_validation.csv`.*

The held-out C is higher than the training C. With 37 deaths in the test set, this difference is within the
range of random variation, so the two values should not be compared closely.

## 5. Discussion

**Interpretation.** Liver function markers drive prognosis in PBC. Higher bilirubin and lower albumin were the
strongest and most reliable signals, consistent with clinical practice. Histologic stage separated survival
clearly on its own, but adding the other variables left stage without a clear independent effect. This is
expected, because stage overlaps with bilirubin and albumin. D-penicillamine showed no benefit. The
confidence interval (0.62 to 1.29) includes both a meaningful benefit and a meaningful harm, so the trial
data cannot rule out either.

**Relation to the real-world problem.** A prognostic model built from routine tests could help clinicians
identify patients who need closer monitoring or earlier transplant referral. The bilirubin-only model achieved
C = 0.82, which is useful. Adding the other variables raised it to 0.86. Any clinical tool would need
validation in current patients, who are treated with ursodeoxycholic acid, a drug that did not exist in this
trial's period. Treatment differences in this dataset therefore reflect a different era of care.

**Limitations.**
1. *Historical data.* The trial ran in 1974 to 1984. Current treatment, including ursodeoxycholic acid, was not
   available then, so survival and risk factors may differ today.
2. *Transplant as censoring.* Treating transplant as censoring assumes transplanted patients would have had the
   same risk as the rest of the cohort. This is probably not true, and it can bias estimates.
3. *Non-randomised subset.* The 106 non-randomised patients have fewer measurements. The main model uses only the
   randomised patients, so it cannot be extended to all 418 without imputation.
4. *Protime estimate.* The protime coefficient is very large and imprecise, and its proportional hazards assumption
   is violated. It should not be used as a clinical threshold.
5. *Small validation set.* The held-out set has only 37 deaths. Concordance estimates have wide uncertainty, and a
   single split is not enough to rank models with confidence.
6. *Single dataset.* The data come from one clinic. External validation in another cohort is needed before use.
7. *Missing stage.* Six patients lack stage and are excluded from stage analyses.

**Real-world implications.** The results support using routine liver function measures in prognosis. They do not
show that D-penicillamine is useless for patients today, since the drug was tested under different conditions and
has since been superseded. Any prognostic tool should be tested on current, prospectively collected data.

## 6. Conclusion and future work

In the Mayo Clinic PBC data, serum bilirubin, albumin, oedema and age were associated with death, and histologic
stage predicted survival on its own. D-penicillamine showed no clear survival benefit. A multivariable model
predicted death in held-out patients with C = 0.86.

Future work should:
- repeat the validation across many random splits and report confidence intervals for C;
- handle transplant as a competing risk rather than as censoring (for example with a Fine-Gray model);
- model the time-varying effects of oedema and protime, for example with a time-dependent coefficient;
- compare with a modern cohort treated with ursodeoxycholic acid;
- use the longitudinal `pbcseq` dataset in the same package to model changes in bilirubin over time.

## 7. Reproducibility

**Requirements.** R 4.3.3 or later. Packages: `survival` (3.5-8 or later), `tidyverse` (2.0.0 or later).

**Steps.** From the project root:

```bash
Rscript code/00_export_data.R   # writes data/pbc.csv from the survival package
Rscript code/01_analysis.R      # analysis, tables and figures
```

The analysis sets `set.seed(2026)` and uses relative paths. It writes tables to `results/`, figures to
`results/figures/`, plus `results/run_log.txt` and `results/session_info.txt`.

**Rendered output.** Plain R scripts; no HTML or PDF render is included.

## 8. References

Cox, D. R. (1972). Regression models and life-tables. *Journal of the Royal Statistical Society: Series B*,
34(2), 187-220.

Fleming, T. R., & Harrington, D. P. (1991). *Counting processes and survival analysis*. Wiley.

Harrell, F. E., Lee, K. L., & Mark, D. B. (1996). Multivariable prognostic models: issues in developing models,
evaluating assumptions and adequacy, and measuring and reducing errors. *Statistics in Medicine*, 15(4), 361-387.

Kaplan, E. L., & Meier, P. (1958). Nonparametric estimation from incomplete observations. *Journal of the American
Statistical Association*, 53(282), 457-481.

Schoenfeld, D. (1982). Partial residuals for the proportional hazards regression model. *Biometrika*, 69(1), 239-241.

Therneau, T. (2024). *A package for survival analysis in R* (R package version 3.5-8).
https://CRAN.R-project.org/package=survival

Therneau, T. M., & Grambsch, P. M. (2000). *Modeling survival data: Extending the Cox model*. Springer.

Wickham, H. (2016). *ggplot2: Elegant graphics for data analysis*. Springer.

Wickham, H., et al. (2019). Welcome to the tidyverse. *Journal of Open Source Software*, 4(43), 1686.
