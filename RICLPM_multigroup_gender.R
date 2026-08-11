library(lavaan)
library(readxl)
library(openxlsx)

data_file <- "C:/R/JYA/data.xlsx"
output_file <- "C:/R/JYA/RI_CLPM_multigroup_gender_12crosspaths_results.xlsx"

# x1-x3 = positive social media emotional responses
# m1-m3 = negative social media emotional responses
# y1-y3 = subjective well-being
# Note: m is a variable name only; it is not treated as a mediator.

required_vars <- c(
  "gender", "age", "ses",
  "x1", "x2", "x3",
  "m1", "m2", "m3",
  "y1", "y2", "y3"
)

dat <- as.data.frame(read_excel(data_file))

missing_vars <- setdiff(required_vars, names(dat))
if (length(missing_vars) > 0) {
  stop("Missing required variables: ", paste(missing_vars, collapse = ", "))
}

analysis_vars <- setdiff(required_vars, "gender")
non_numeric_vars <- analysis_vars[
  !vapply(dat[analysis_vars], is.numeric, logical(1))
]
if (length(non_numeric_vars) > 0) {
  stop(
    "These analysis variables must be numeric: ",
    paste(non_numeric_vars, collapse = ", ")
  )
}

# The grouping variable itself cannot be handled by FIML.
gender_text <- trimws(as.character(dat$gender))
missing_gender <- is.na(dat$gender) | gender_text == ""
n_missing_gender <- sum(missing_gender)

if (n_missing_gender > 0) {
  message("Excluded ", n_missing_gender, " cases with missing gender.")
  dat <- dat[!missing_gender, , drop = FALSE]
  gender_text <- gender_text[!missing_gender]
}

dat$gender_group <- factor(gender_text)
group_levels <- levels(dat$gender_group)

if (length(group_levels) != 2) {
  stop(
    "gender must contain exactly two non-missing groups; found: ",
    paste(group_levels, collapse = ", ")
  )
}

group_counts <- table(dat$gender_group)
group_info <- data.frame(
  lavaan_group = 1:2,
  original_gender_value = group_levels,
  n = as.integer(group_counts[group_levels]),
  excluded_missing_gender_total = n_missing_gender,
  stringsAsFactors = FALSE
)

############################################################
# Common RI-CLPM structure: nuisance parameters are free    #
# across gender; all within-person paths are time-specific. #
############################################################

model_before_paths <- '
  RIx =~ 1*x1 + 1*x2 + 1*x3
  RIm =~ 1*m1 + 1*m2 + 1*m3
  RIy =~ 1*y1 + 1*y2 + 1*y3

  RIx ~ c(rix_age_g1, rix_age_g2)*age + c(rix_ses_g1, rix_ses_g2)*ses
  RIm ~ c(rim_age_g1, rim_age_g2)*age + c(rim_ses_g1, rim_ses_g2)*ses
  RIy ~ c(riy_age_g1, riy_age_g2)*age + c(riy_ses_g1, riy_ses_g2)*ses

  x1 ~ c(mux1_g1, mux1_g2)*1
  x2 ~ c(mux2_g1, mux2_g2)*1
  x3 ~ c(mux3_g1, mux3_g2)*1
  m1 ~ c(mum1_g1, mum1_g2)*1
  m2 ~ c(mum2_g1, mum2_g2)*1
  m3 ~ c(mum3_g1, mum3_g2)*1
  y1 ~ c(muy1_g1, muy1_g2)*1
  y2 ~ c(muy2_g1, muy2_g2)*1
  y3 ~ c(muy3_g1, muy3_g2)*1

  RIx ~~ c(vrix_g1, vrix_g2)*RIx
  RIm ~~ c(vrim_g1, vrim_g2)*RIm
  RIy ~~ c(vriy_g1, vriy_g2)*RIy
  RIx ~~ c(crixm_g1, crixm_g2)*RIm
  RIx ~~ c(crixy_g1, crixy_g2)*RIy
  RIm ~~ c(crimy_g1, crimy_g2)*RIy

  wx1 =~ 1*x1
  wx2 =~ 1*x2
  wx3 =~ 1*x3
  wm1 =~ 1*m1
  wm2 =~ 1*m2
  wm3 =~ 1*m3
  wy1 =~ 1*y1
  wy2 =~ 1*y2
  wy3 =~ 1*y3
'

# Free model: all 18 within-person paths may differ by gender.
free_paths <- '
  wx2 ~ c(ax1_g1, ax1_g2)*wx1 + c(mx1_g1, mx1_g2)*wm1 + c(yx1_g1, yx1_g2)*wy1
  wm2 ~ c(xm1_g1, xm1_g2)*wx1 + c(am1_g1, am1_g2)*wm1 + c(ym1_g1, ym1_g2)*wy1
  wy2 ~ c(xy1_g1, xy1_g2)*wx1 + c(my1_g1, my1_g2)*wm1 + c(ay1_g1, ay1_g2)*wy1

  wx3 ~ c(ax2_g1, ax2_g2)*wx2 + c(mx2_g1, mx2_g2)*wm2 + c(yx2_g1, yx2_g2)*wy2
  wm3 ~ c(xm2_g1, xm2_g2)*wx2 + c(am2_g1, am2_g2)*wm2 + c(ym2_g1, ym2_g2)*wy2
  wy3 ~ c(xy2_g1, xy2_g2)*wx2 + c(my2_g1, my2_g2)*wm2 + c(ay2_g1, ay2_g2)*wy2
'

# Constrained model: only the 12 cross-lagged paths are equal across gender.
# The 6 autoregressive paths remain freely estimated in both groups.
equal_paths <- '
  wx2 ~ c(ax1_g1, ax1_g2)*wx1 + c(mx1, mx1)*wm1 + c(yx1, yx1)*wy1
  wm2 ~ c(xm1, xm1)*wx1 + c(am1_g1, am1_g2)*wm1 + c(ym1, ym1)*wy1
  wy2 ~ c(xy1, xy1)*wx1 + c(my1, my1)*wm1 + c(ay1_g1, ay1_g2)*wy1

  wx3 ~ c(ax2_g1, ax2_g2)*wx2 + c(mx2, mx2)*wm2 + c(yx2, yx2)*wy2
  wm3 ~ c(xm2, xm2)*wx2 + c(am2_g1, am2_g2)*wm2 + c(ym2, ym2)*wy2
  wy3 ~ c(xy2, xy2)*wx2 + c(my2, my2)*wm2 + c(ay2_g1, ay2_g2)*wy2
'

model_after_paths <- '
  wx1 ~~ c(vwx1_g1, vwx1_g2)*wx1
  wx2 ~~ c(vx_g1, vx_g2)*wx2
  wx3 ~~ c(vx_g1, vx_g2)*wx3
  wm1 ~~ c(vwm1_g1, vwm1_g2)*wm1
  wm2 ~~ c(vm_g1, vm_g2)*wm2
  wm3 ~~ c(vm_g1, vm_g2)*wm3
  wy1 ~~ c(vwy1_g1, vwy1_g2)*wy1
  wy2 ~~ c(vy_g1, vy_g2)*wy2
  wy3 ~~ c(vy_g1, vy_g2)*wy3

  wx1 ~~ c(cxm1_g1, cxm1_g2)*wm1
  wx1 ~~ c(cxy1_g1, cxy1_g2)*wy1
  wm1 ~~ c(cmy1_g1, cmy1_g2)*wy1
  wx2 ~~ c(cxm2_g1, cxm2_g2)*wm2
  wx2 ~~ c(cxy2_g1, cxy2_g2)*wy2
  wm2 ~~ c(cmy2_g1, cmy2_g2)*wy2
  wx3 ~~ c(cxm3_g1, cxm3_g2)*wm3
  wx3 ~~ c(cxy3_g1, cxy3_g2)*wy3
  wm3 ~~ c(cmy3_g1, cmy3_g2)*wy3
'

model_free <- paste(model_before_paths, free_paths, model_after_paths, sep = "\n")
model_equal <- paste(model_before_paths, equal_paths, model_after_paths, sep = "\n")

fit_riclpm <- function(model) {
  lavaan(
    model = model,
    data = dat,
    group = "gender_group",
    group.label = group_levels,
    estimator = "MLR",
    missing = "FIML",
    fixed.x = TRUE,
    int.ov.free = FALSE,
    int.lv.free = FALSE,
    auto.fix.first = FALSE,
    auto.fix.single = FALSE,
    auto.cov.lv.x = FALSE,
    auto.cov.y = FALSE,
    auto.var = FALSE
  )
}

fit_free <- fit_riclpm(model_free)
fit_equal <- fit_riclpm(model_equal)

check_fit <- function(fit, model_name) {
  if (!isTRUE(lavInspect(fit, "converged"))) {
    stop(model_name, " did not converge.")
  }
  if (!isTRUE(lavInspect(fit, "post.check"))) {
    stop(model_name, " has an inadmissible solution (for example, a negative variance).")
  }
}

check_fit(fit_free, "Gender-free model")
check_fit(fit_equal, "Cross-lag-constrained model")

############################################################
# Model fit and robust likelihood-ratio comparison          #
############################################################

fit_measures <- c(
  "chisq.scaled", "df.scaled", "pvalue.scaled",
  "cfi.robust", "rmsea.robust", "srmr", "aic", "bic"
)

fit_row <- function(fit, model_name) {
  values <- fitMeasures(fit, fit_measures)
  out <- as.data.frame(as.list(values), check.names = FALSE)
  cbind(Model = model_name, out)
}

fit_comparison <- rbind(
  fit_row(fit_free, "Gender-free"),
  fit_row(fit_equal, "Cross-lag-constrained")
)
rownames(fit_comparison) <- NULL

df_difference <- unname(
  fitMeasures(fit_equal, "df") - fitMeasures(fit_free, "df")
)
if (!isTRUE(all.equal(df_difference, 12))) {
  stop("Expected a 12-df difference, but found ", df_difference, ".")
}

robust_lrt <- as.data.frame(
  lavTestLRT(fit_free, fit_equal, method = "default"),
  check.names = FALSE
)
robust_lrt <- cbind(Model = rownames(robust_lrt), robust_lrt)
rownames(robust_lrt) <- NULL

p_columns <- grep("^Pr\\(", names(robust_lrt), value = TRUE)
if (length(p_columns) == 0) {
  stop("Could not find the robust likelihood-ratio p-value.")
}
global_p <- robust_lrt[[p_columns[length(p_columns)]]][nrow(robust_lrt)]
global_significant <- !is.na(global_p) && global_p < 0.05

############################################################
# Cross-lagged path Wald tests with BH-FDR correction       #
############################################################

path_map <- data.frame(
  path_id = c(
    "m1_to_x2", "y1_to_x2", "x1_to_m2", "y1_to_m2", "x1_to_y2", "m1_to_y2",
    "m2_to_x3", "y2_to_x3", "x2_to_m3", "y2_to_m3", "x2_to_y3", "m2_to_y3"
  ),
  path = c(
    "wm1 -> wx2", "wy1 -> wx2", "wx1 -> wm2", "wy1 -> wm2", "wx1 -> wy2", "wm1 -> wy2",
    "wm2 -> wx3", "wy2 -> wx3", "wx2 -> wm3", "wy2 -> wm3", "wx2 -> wy3", "wm2 -> wy3"
  ),
  label_g1 = c(
    "mx1_g1", "yx1_g1", "xm1_g1", "ym1_g1", "xy1_g1", "my1_g1",
    "mx2_g1", "yx2_g1", "xm2_g1", "ym2_g1", "xy2_g1", "my2_g1"
  ),
  label_g2 = c(
    "mx1_g2", "yx1_g2", "xm1_g2", "ym1_g2", "xy1_g2", "my1_g2",
    "mx2_g2", "yx2_g2", "xm2_g2", "ym2_g2", "xy2_g2", "my2_g2"
  ),
  stringsAsFactors = FALSE
)

pe_free <- parameterEstimates(fit_free, ci = TRUE, standardized = FALSE)

get_estimate <- function(label) {
  rows <- pe_free[pe_free$label == label & pe_free$op == "~", , drop = FALSE]
  if (nrow(rows) == 0) {
    stop("Parameter label not found: ", label)
  }
  c(estimate = rows$est[1], se = rows$se[1])
}

wald_rows <- lapply(seq_len(nrow(path_map)), function(i) {
  est_g1 <- get_estimate(path_map$label_g1[i])
  est_g2 <- get_estimate(path_map$label_g2[i])
  wald <- lavTestWald(
    fit_free,
    constraints = paste(path_map$label_g1[i], "==", path_map$label_g2[i])
  )

  data.frame(
    path_id = path_map$path_id[i],
    path = path_map$path[i],
    group_1_value = group_levels[1],
    estimate_group_1 = unname(est_g1["estimate"]),
    se_group_1 = unname(est_g1["se"]),
    group_2_value = group_levels[2],
    estimate_group_2 = unname(est_g2["estimate"]),
    se_group_2 = unname(est_g2["se"]),
    Wald_W = unname(wald$stat),
    df = unname(wald$df),
    p_value = unname(wald$p.value),
    stringsAsFactors = FALSE
  )
})

wald_fdr <- do.call(rbind, wald_rows)
rownames(wald_fdr) <- NULL

if (nrow(wald_fdr) != 12 || anyDuplicated(wald_fdr$path_id)) {
  stop("The Wald table must contain exactly 12 unique cross-lagged paths.")
}
if (anyNA(wald_fdr$p_value)) {
  stop("At least one Wald test returned a missing p-value.")
}

wald_fdr$p_fdr_bh <- p.adjust(wald_fdr$p_value, method = "BH")
wald_fdr$global_LRT_p <- global_p
wald_fdr$global_LRT_significant_05 <- global_significant
wald_fdr$gender_difference_supported <-
  global_significant & wald_fdr$p_fdr_bh < 0.05

############################################################
# Export                                                     #
############################################################

write.xlsx(
  list(
    group_info = group_info,
    fit_comparison = fit_comparison,
    robust_LRT = robust_lrt,
    crosslag_Wald_FDR = wald_fdr
  ),
  file = output_file,
  overwrite = TRUE
)

cat("\nGroup information:\n")
print(group_info)
cat("\nModel fit comparison:\n")
print(fit_comparison)
cat("\nRobust likelihood-ratio test:\n")
print(robust_lrt)
cat("\nCross-lagged path Wald tests with BH-FDR:\n")
print(wald_fdr)
cat("\nResults saved to: ", output_file, "\n", sep = "")

openxlsx::write.xlsx( wald_fdr,  file = "C:/R/JYA/RI_CLPM_gender_robust_LRT_敏感性分析.xlsx",
  rowNames = FALSE,
  overwrite = TRUE
)
