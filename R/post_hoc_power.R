# Non-inferiority analysis of seroconversion: for each state, age group and
# poliovirus type, estimate the paired risk difference in neutralising-antibody
# positivity between baseline and follow-up, with a skewness-corrected score
# confidence interval and a one-sided non-inferiority test against two margins
# (0.10 and 0.30). Output: output/non_inferiority_analysis.docx

source(here::here("R", "libraries.R"))

load(here(path("data","anon_neut.Rdata")))
load(here(path("data","anon_vacc.Rdata")))


# Paired 2x2 counts per stratum (positive = neut_pos), baseline (b) by follow-up (f):
#   a = positive at both, b = negative at baseline & positive at follow-up,
#   c = positive at baseline & negative at follow-up, d = negative at both
data<-
df_neut|>
  pivot_wider(names_from=time,values_from=c(neut_pos,neut_pos_SA,neut),id_cols=c(pv,id,sample_print))|>
  left_join(df_vacc|>select(age_cat,state_print,nopv2_sia,bopv_ri_3plus_num,ipv_1plus_num,ri_info_source,id),by="id")|>
  summarise(.by=c(state_print,age_cat,pv),
a=sum(neut_pos_b&neut_pos_f,na.rm=T),
b=sum(!neut_pos_b&neut_pos_f,na.rm=T),
c=sum(neut_pos_b&!neut_pos_f,na.rm=T),
d=sum(!neut_pos_b&!neut_pos_f,na.rm=T))


# install.packages("ratesci")

# Candidate non-inferiority margins (risk difference)
margins <- c(0.10, 0.30)

# Cross each stratum with each candidate margin
data_by_margin <- tidyr::crossing(data, margin = margins)

# For each stratum and margin, fit the paired risk-difference score test and
# collect the estimate, 95% CI and p-value
results <- purrr::pmap(data_by_margin, function(state_print, age_cat, pv, a, b, c, d, margin, ...) {
  x <- c(a = a, b = b, c = c, d = d)

  fit <- scorepairci(
    x = x,
    contrast = "RD",
    level = 0.95,
    theta0 = margin
  )

  data.frame(
    state_print  = state_print,
    age_cat      = age_cat,
    pv           = pv,
    margin       = margin,
    RD           = fit$estimates[, "est"],
    CI_lower     = fit$estimates[, "lower"],
    CI_upper     = fit$estimates[, "upper"],
    # One-sided NI p-value: H0: RD >= margin vs H1: RD < margin
    p_value_NI   = fit$pval[, "pval_left"]
  )
}) |>
  dplyr::bind_rows() |>
  # Bonferroni correction for multiple testing across the 12 strata,
  # applied separately within each margin
  mutate(p_value_NI_adj = p.adjust(p_value_NI, method = "bonferroni"), .by = margin)

# One row per stratum, with p-values for each margin in separate columns
results_wide <-
  results |>
  pivot_wider(
    id_cols = c(state_print, age_cat, pv, RD, CI_lower, CI_upper),
    names_from = margin,
    values_from = c(p_value_NI, p_value_NI_adj)
  )

# Formatted results table, grouped by poliovirus type
ni_table <-
  results_wide |>
  mutate(pv = paste0("Poliovirus type ", pv)) |>
  arrange(pv, state_print, age_cat) |>
  gt(groupname_col = "pv") |>
  cols_label(
    state_print         = "State",
    age_cat             = "Age group",
    RD                  = "Risk difference",
    CI_lower            = "Lower",
    CI_upper            = "Upper",
    p_value_NI_0.1      = md("*p*-value"),
    p_value_NI_adj_0.1  = md("Adjusted (Bonferroni)"),
    p_value_NI_0.3      = md("*p*-value"),
    p_value_NI_adj_0.3  = md("Adjusted (Bonferroni)")
  ) |>
  tab_spanner(
    label   = "95% CI",
    columns = c(CI_lower, CI_upper)
  ) |>
  tab_spanner(
    label   = "NI margin = 0.10",
    columns = c(p_value_NI_0.1, p_value_NI_adj_0.1)
  ) |>
  tab_spanner(
    label   = "NI margin = 0.30",
    columns = c(p_value_NI_0.3, p_value_NI_adj_0.3)
  ) |>
  tab_header(
    title    = "Non-inferiority analysis: seroconversion risk difference (follow-up to baseline)",
    subtitle = "Skewness-corrected score method; margins = 0.10 and 0.30"
  ) |>
  fmt_number(columns = c(RD, CI_lower, CI_upper), decimals = 2) |>
  # Show p-values rounding to 0.000 as "<0.001" instead
  fmt(
    columns = c(p_value_NI_0.1, p_value_NI_adj_0.1, p_value_NI_0.3, p_value_NI_adj_0.3),
    fns = function(x) ifelse(x < 0.0005, "<0.001", formatC(x, format = "f", digits = 3))
  ) |>
  # Bold Bonferroni-adjusted p-values that fail to reach significance (> 0.05)
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = p_value_NI_adj_0.1,
      rows = p_value_NI_adj_0.1 > 0.05
    )
  ) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = p_value_NI_adj_0.3,
      rows = p_value_NI_adj_0.3 > 0.05
    )
  )

ni_table

gtsave(ni_table, filename = here(path("output", "non_inferiority_analysis.docx")))

