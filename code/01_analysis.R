# =============================================================================
# Survival in primary biliary cholangitis (PBC): which baseline measurements
# predict death, and does D-penicillamine change survival?
#
# Data: Mayo Clinic PBC trial, 1974-1984 (418 patients), from the 'survival'
# package, exported to data/pbc.csv by code/00_export_data.R.
#
# Outcome: death. Liver transplantation (status 1) is treated as censoring at
# transplant. Patients alive at the end of follow-up are censored.
#
# Run from the project root:  Rscript code/01_analysis.R
# =============================================================================

suppressPackageStartupMessages({
  library(survival)   
  library(tidyverse)  
})

set.seed(2026)
FIG <- "C:/Users/User/Desktop/Gawesha/pbc-liver-survival-r/pbc-liver-survival-r/results/figures"
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
theme_set(theme_minimal(base_size = 11))
pal <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7")

# ---- 1. Import ------------------------------------------------------------------
raw <- read_csv("C:/Users/User/Desktop/Gawesha/pbc-liver-survival-r/pbc-liver-survival-r/data/pbc.csv", show_col_types = FALSE)

# ---- 2. Quality control ---------------------------------------------------------
qc <- tibble(
  n_patients      = nrow(raw),
  n_randomised    = sum(!is.na(raw$trt)),
  n_not_randomised = sum(is.na(raw$trt)),
  n_deaths        = sum(raw$status == 2),
  n_transplants   = sum(raw$status == 1),
  n_censored      = sum(raw$status == 0),
  median_followup_years = round(median(raw$time) / 365.25, 2)
)
write_csv(qc, "results/qc_summary.csv")
print(qc)

missing <- tibble(variable = names(raw),
                  n_missing = colSums(is.na(raw)),
                  pct_missing = round(100 * colMeans(is.na(raw)), 1)) |>
  filter(n_missing > 0) |> arrange(desc(n_missing))
write_csv(missing, "results/missingness.csv")

# ---- 3. Cleaning and derived variables ----------------------------------------------
dat <- raw |>
  mutate(
    years    = time / 365.25,
    event    = as.integer(status == 2),     # 1 = died, 0 = alive or transplanted
    trt      = factor(trt, levels = c(2, 1),
                      labels = c("Placebo", "D-penicillamine")),
    sex      = factor(sex, levels = c("f", "m"), labels = c("Female", "Male")),
    stage    = factor(stage, levels = 1:4, labels = paste("Stage", 1:4)),
    ascites  = factor(ascites, levels = c(0, 1), labels = c("No", "Yes")),
    hepato   = factor(hepato, levels = c(0, 1), labels = c("No", "Yes")),
    spiders  = factor(spiders, levels = c(0, 1), labels = c("No", "Yes")),
    log_bili = log(bili),         # bilirubin is strongly right-skewed
    log_protime = log(protime)
  )

# The randomised trial subset (312 patients) has the full set of clinical measures.
# The 106 non-randomised patients have fewer measures, so the main model uses the
# randomised subset only.
trial <- dat |> filter(!is.na(trt))

model_vars <- c("trt", "age", "sex", "ascites", "hepato", "spiders", "edema",
                "log_bili", "albumin", "log_protime", "stage")
model_dat <- trial |> drop_na(all_of(model_vars))
cat(sprintf("Randomised patients: %d; complete for model: %d (dropped %d)\n",
            nrow(trial), nrow(model_dat), nrow(trial) - nrow(model_dat)))
write_csv(tibble(analysis = c("randomised trial", "model (complete cases)"),
                 n = c(nrow(trial), nrow(model_dat)),
                 deaths = c(sum(trial$event), sum(model_dat$event))),
          "results/analysis_sample.csv")

# Baseline characteristics by treatment (randomised subset)
baseline <- trial |>
  group_by(trt) |>
  summarise(n = n(),
            deaths = sum(event),
            age_mean = round(mean(age), 1),
            female_pct = round(100 * mean(sex == "Female"), 1),
            bili_median = round(median(bili), 2),
            albumin_median = round(median(albumin), 2),
            protime_median = round(median(protime), 2),
            .groups = "drop")
write_csv(baseline, "results/baseline_by_treatment.csv")
print(baseline)

# ---- 4. Exploratory analysis ---------------------------------------------------------
# Fig 1: follow-up time by outcome
p1 <- dat |>
  mutate(outcome = case_when(status == 2 ~ "Died",
                             status == 1 ~ "Transplanted (censored)",
                             TRUE ~ "Alive (censored)")) |>
  ggplot(aes(years, fill = outcome)) +
  geom_histogram(bins = 30, colour = "white", position = "stack") +
  scale_fill_manual(values = c("Died" = pal[2],
                               "Transplanted (censored)" = pal[1],
                               "Alive (censored)" = "grey70")) +
  labs(title = "Fig 1. Follow-up time by outcome",
       x = "Years from registration", y = "Number of patients", fill = NULL)
ggsave(file.path(FIG, "fig1_followup.png"), p1, width = 7, height = 4.5, dpi = 150)

# Fig 2: bilirubin by outcome (log scale)
p2 <- ggplot(dat, aes(factor(event, labels = c("Alive or transplanted", "Died")),
                      bili, fill = factor(event))) +
  geom_boxplot(outlier.alpha = 0.4) +
  scale_y_log10() +
  scale_fill_manual(values = c("0" = "grey70", "1" = pal[2]), guide = "none") +
  labs(title = "Fig 2. Serum bilirubin by outcome",
       x = NULL, y = "Bilirubin (mg/dl, log scale)")
ggsave(file.path(FIG, "fig2_bilirubin_by_outcome.png"), p2, width = 6, height = 4.5, dpi = 150)

# ---- 5. Kaplan-Meier analysis ---------------------------------------------------------
km_tbl <- function(fit) {
  s <- summary(fit)
  tibble(time = s$time, surv = s$surv, lower = s$lower, upper = s$upper,
         group = if (is.null(s$strata)) "All" else sub("^.*=", "", as.character(s$strata)))
}

km_all <- survfit(Surv(years, event) ~ 1, data = dat)
km_median <- summary(km_all)$table["median"]
cat(sprintf("Overall median survival: %.2f years\n", km_median))

# By histologic stage (all patients with stage recorded)
km_stage_dat <- dat |> filter(!is.na(stage))
km_stage <- survfit(Surv(years, event) ~ stage, data = km_stage_dat)
lr_stage <- survdiff(Surv(years, event) ~ stage, data = km_stage_dat)
lr_stage_p <- 1 - pchisq(lr_stage$chisq, df = length(lr_stage$n) - 1)
cat(sprintf("Log-rank test, stage: p = %.3g\n", lr_stage_p))

p3 <- km_tbl(km_stage) |>
  ggplot(aes(time, surv, colour = group, fill = group)) +
  geom_step(linewidth = 0.8) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1, colour = NA) +
  scale_colour_manual(values = pal) + scale_fill_manual(values = pal) +
  labs(title = "Fig 3. Survival by histologic stage",
       subtitle = sprintf("Kaplan-Meier with 95%% CI; log-rank p = %.3g", lr_stage_p),
       x = "Years from registration", y = "Probability of survival", colour = NULL, fill = NULL) +
  ylim(0, 1) + theme(legend.position = "bottom")
ggsave(file.path(FIG, "fig3_km_stage.png"), p3, width = 7, height = 5, dpi = 150)

# By treatment (randomised trial only)
km_trt <- survfit(Surv(years, event) ~ trt, data = trial)
lr_trt <- survdiff(Surv(years, event) ~ trt, data = trial)
lr_trt_p <- 1 - pchisq(lr_trt$chisq, df = 1)
cat(sprintf("Log-rank test, treatment: p = %.3g\n", lr_trt_p))

p4 <- km_tbl(km_trt) |>
  ggplot(aes(time, surv, colour = group, fill = group)) +
  geom_step(linewidth = 0.8) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1, colour = NA) +
  scale_colour_manual(values = c(pal[1], pal[2])) +
  scale_fill_manual(values = c(pal[1], pal[2])) +
  labs(title = "Fig 4. Survival by randomised treatment",
       subtitle = sprintf("Randomised trial subset (n = %d); log-rank p = %.2f",
                          nrow(trial), lr_trt_p),
       x = "Years from registration", y = "Probability of survival", colour = NULL, fill = NULL) +
  ylim(0, 1) + theme(legend.position = "bottom")
ggsave(file.path(FIG, "fig4_km_treatment.png"), p4, width = 7, height = 5, dpi = 150)

write_csv(tibble(comparison = c("stage", "treatment"),
                 test = "log-rank",
                 p_value = c(lr_stage_p, lr_trt_p)),
          "results/logrank_tests.csv")

# ---- 6. Cox proportional hazards models --------------------------------------------
# Univariable: each variable alone (model sample)
uni <- map_dfr(model_vars, function(v) {
  f <- as.formula(paste("Surv(years, event) ~", v))
  fit <- coxph(f, data = model_dat)
  ci <- summary(fit)$conf.int
  tibble(term = rownames(ci)[1], HR = ci[1, 1], lower = ci[1, 3], upper = ci[1, 4],
         p_value = summary(fit)$coefficients[1, 5])
})
write_csv(uni |> mutate(across(where(is.numeric), ~ signif(.x, 3))),
          "results/cox_univariable.csv")

# Multivariable: all candidate variables together
cox_multi <- coxph(Surv(years, event) ~ trt + age + sex + ascites + hepato + spiders +
                     edema + log_bili + albumin + log_protime + stage,
                   data = model_dat)
print(summary(cox_multi))

multi_tbl <- tibble(
  term = rownames(summary(cox_multi)$conf.int),
  HR = summary(cox_multi)$conf.int[, 1],
  lower = summary(cox_multi)$conf.int[, 3],
  upper = summary(cox_multi)$conf.int[, 4],
  p_value = summary(cox_multi)$coefficients[, 5]
) |> mutate(across(where(is.numeric), ~ signif(.x, 3)))
write_csv(multi_tbl, "results/cox_multivariable.csv")

# Fig 5: forest plot of adjusted hazard ratios
p5 <- multi_tbl |>
  mutate(term = fct_reorder(term, HR)) |>
  ggplot(aes(HR, term)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, colour = "grey40") +
  geom_point(colour = pal[1], size = 2) +
  scale_x_log10() +
  labs(title = "Fig 5. Adjusted hazard ratios for death",
       subtitle = "Multivariable Cox model, randomised subset (complete cases); log scale",
       x = "Hazard ratio (95% CI)", y = NULL)
ggsave(file.path(FIG, "fig5_forest.png"), p5, width = 7, height = 5.5, dpi = 150)

# ---- 7. Proportional hazards assumption ----------------------------------------------
ph <- cox.zph(cox_multi)
ph_tbl <- tibble(term = rownames(ph$table), chisq = ph$table[, 1],
                 df = ph$table[, 2], p_value = ph$table[, 3]) |>
  mutate(across(where(is.numeric), ~ signif(.x, 3)))
write_csv(ph_tbl, "results/ph_test.csv")
print(ph_tbl)

png(file.path(FIG, "fig6_schoenfeld_log_bili.png"), width = 1000, height = 700, res = 150)
plot(ph[which(rownames(ph$table) == "log_bili")],
     main = "Fig 6. Schoenfeld residuals: log bilirubin", resid = FALSE)
abline(h = 0, lty = 2)
invisible(dev.off())

# ---- 8. Validation: 70/30 split with stratification by outcome --------------------------
# The model is fitted on the training set only, then scored on the held-out set
# with the concordance index (C), where 0.5 is chance and 1 is perfect ranking.
train_idx <- unlist(lapply(split(seq_len(nrow(model_dat)), model_dat$event),
                           function(i) sample(i, size = round(0.7 * length(i)))))
train <- model_dat[train_idx, ]
test  <- model_dat[-train_idx, ]

fit_full_train <- coxph(Surv(years, event) ~ trt + age + sex + ascites + hepato + spiders +
                          edema + log_bili + albumin + log_protime + stage, data = train)
fit_bili_train <- coxph(Surv(years, event) ~ log_bili, data = train)

c_index <- function(fit, newdata) {
  lp <- predict(fit, newdata = newdata, type = "lp")
  # Higher linear predictor = higher risk = shorter survival, so reverse = TRUE
  concordance(Surv(years, event) ~ lp, data = newdata, reverse = TRUE)$concordance
}

validation <- tibble(
  model = c("Multivariable Cox (all candidate variables)", "Bilirubin only"),
  C_train = c(concordance(fit_full_train)$concordance, concordance(fit_bili_train)$concordance),
  C_test  = c(c_index(fit_full_train, test), c_index(fit_bili_train, test)),
  n_train = nrow(train), n_test = nrow(test), deaths_test = sum(test$event)
) |> mutate(across(starts_with("C_"), ~ round(as.numeric(.x), 3)))
write_csv(validation, "results/model_validation.csv")
print(validation)

# ---- 9. Session information ----------------------------------------------------------
writeLines(capture.output(sessionInfo()), "results/session_info.txt")
cat("Analysis complete. Outputs written to results/.\n")
