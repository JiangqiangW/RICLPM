library(lavaan)
library(readxl)
library(openxlsx)

setwd("C:/R/JYA")

# Required columns:
# gender, age, ses, x1, x2, x3, m1, m2, m3, y1, y2, y3
# x1-x3 = positive social media emotional responses
# m1-m3 = negative social media emotional responses
# y1-y3 = subjective well-being
# Note: m is a variable name only; it is not treated as a mediator.

dat <- read_excel("JYA_C01.xlsx")
dat <- as.data.frame(dat)

lav_options <- list(
  data = dat,
  estimator = "MLR",
  missing = "FIML",
  int.ov.free = FALSE,
  int.lv.free = FALSE,
  auto.fix.first = FALSE,
  auto.fix.single = FALSE,
  auto.cov.lv.x = FALSE,
  auto.cov.y = FALSE,
  auto.var = FALSE)

fit_riclpm <- function(model) {
  do.call(lavaan, c(list(model = model), lav_options))
}

fit_table <- function(fits) {
  measures <- c(
    "chisq.scaled", "df.scaled", "pvalue.scaled",
    "cfi.robust", "aic", "bic",
    "rmsea.robust", "rmsea.ci.lower.robust", "rmsea.ci.upper.robust",
    "srmr"
  )

  out <- do.call(rbind, lapply(fits, function(fit) fitMeasures(fit, measures)))
  out <- as.data.frame(out)
  out <- cbind(Model = names(fits), out)

  names(out) <- c(
    "Model", "ChiSquare", "DF", "PValue", "CFI", "AIC", "BIC",
    "RMSEA", "RMSEA_Lower", "RMSEA_Upper", "SRMR"
  )
  rownames(out) <- NULL
  out
}

##########################################################
### Model 1: no lagged-path constraints                 ###
##########################################################

riclpm_model1 <- '
  RIx =~ 1*x1 + 1*x2 + 1*x3
  RIm =~ 1*m1 + 1*m2 + 1*m3
  RIy =~ 1*y1 + 1*y2 + 1*y3
  

  
  x1 ~ mux1*1
  x2 ~ mux2*1
  x3 ~ mux3*1
  m1 ~ mum1*1
  m2 ~ mum2*1
  m3 ~ mum3*1
  y1 ~ muy1*1
  y2 ~ muy2*1
  y3 ~ muy3*1

  RIx ~~ RIx
  RIm ~~ RIm
  RIy ~~ RIy
  RIx ~~ RIm
  RIx ~~ RIy
  RIm ~~ RIy

  wx1 =~ 1*x1
  wx2 =~ 1*x2
  wx3 =~ 1*x3
  wm1 =~ 1*m1
  wm2 =~ 1*m2
  wm3 =~ 1*m3
  wy1 =~ 1*y1
  wy2 =~ 1*y2
  wy3 =~ 1*y3

  wx2 ~ ax1*wx1 + mx1*wm1 + yx1*wy1
  wm2 ~ xm1*wx1 + am1*wm1 + ym1*wy1
  wy2 ~ xy1*wx1 + my1*wm1 + ay1*wy1

  wx3 ~ ax2*wx2 + mx2*wm2 + yx2*wy2
  wm3 ~ xm2*wx2 + am2*wm2 + ym2*wy2
  wy3 ~ xy2*wx2 + my2*wm2 + ay2*wy2

  wx1 ~~ wx1
  wx2 ~~ vx*wx2
  wx3 ~~ vx*wx3
  wm1 ~~ wm1
  wm2 ~~ vm*wm2
  wm3 ~~ vm*wm3
  wy1 ~~ wy1
  wy2 ~~ vy*wy2
  wy3 ~~ vy*wy3

  wx1 ~~ wm1
  wx1 ~~ wy1
  wm1 ~~ wy1
  wx2 ~~ cxm2*wm2
  wx2 ~~ cxy2*wy2
  wm2 ~~ cmy2*wy2
  wx3 ~~ cxm3*wm3
  wx3 ~~ cxy3*wy3
  wm3 ~~ cmy3*wy3
'

##########################################################
### Model 2: autoregressive paths constrained equal      ###
##########################################################

riclpm_model2 <- '
  RIx =~ 1*x1 + 1*x2 + 1*x3
  RIm =~ 1*m1 + 1*m2 + 1*m3
  RIy =~ 1*y1 + 1*y2 + 1*y3

# RIx ~ gender + age + ses
# RIm ~ gender + age + ses
# RIy ~ gender + age + ses

  
  x1 ~ mux1*1
  x2 ~ mux2*1
  x3 ~ mux3*1
  m1 ~ mum1*1
  m2 ~ mum2*1
  m3 ~ mum3*1
  y1 ~ muy1*1
  y2 ~ muy2*1
  y3 ~ muy3*1

  RIx ~~ RIx
  RIm ~~ RIm
  RIy ~~ RIy
  RIx ~~ RIm
  RIx ~~ RIy
  RIm ~~ RIy

  wx1 =~ 1*x1
  wx2 =~ 1*x2
  wx3 =~ 1*x3
  wm1 =~ 1*m1
  wm2 =~ 1*m2
  wm3 =~ 1*m3
  wy1 =~ 1*y1
  wy2 =~ 1*y2
  wy3 =~ 1*y3

  wx2 ~ ax*wx1 + mx1*wm1 + yx1*wy1
  wm2 ~ xm1*wx1 + am*wm1 + ym1*wy1
  wy2 ~ xy1*wx1 + my1*wm1 + ay*wy1

  wx3 ~ ax*wx2 + mx2*wm2 + yx2*wy2
  wm3 ~ xm2*wx2 + am*wm2 + ym2*wy2
  wy3 ~ xy2*wx2 + my2*wm2 + ay*wy2

  wx1 ~~ wx1
  wx2 ~~ vx*wx2
  wx3 ~~ vx*wx3
  wm1 ~~ wm1
  wm2 ~~ vm*wm2
  wm3 ~~ vm*wm3
  wy1 ~~ wy1
  wy2 ~~ vy*wy2
  wy3 ~~ vy*wy3

  wx1 ~~ wm1
  wx1 ~~ wy1
  wm1 ~~ wy1
  wx2 ~~ cxm2*wm2
  wx2 ~~ cxy2*wy2
  wm2 ~~ cmy2*wy2
  wx3 ~~ cxm3*wm3
  wx3 ~~ cxy3*wy3
  wm3 ~~ cmy3*wy3
'

##########################################################
### Model 3: autoregressive and cross-lagged paths equal ###
##########################################################

riclpm_model3 <- '
  RIx =~ 1*x1 + 1*x2 + 1*x3
  RIm =~ 1*m1 + 1*m2 + 1*m3
  RIy =~ 1*y1 + 1*y2 + 1*y3



  x1 ~ mux1*1
  x2 ~ mux2*1
  x3 ~ mux3*1
  m1 ~ mum1*1
  m2 ~ mum2*1
  m3 ~ mum3*1
  y1 ~ muy1*1
  y2 ~ muy2*1
  y3 ~ muy3*1

  RIx ~~ RIx
  RIm ~~ RIm
  RIy ~~ RIy
  RIx ~~ RIm
  RIx ~~ RIy
  RIm ~~ RIy

  wx1 =~ 1*x1
  wx2 =~ 1*x2
  wx3 =~ 1*x3
  wm1 =~ 1*m1
  wm2 =~ 1*m2
  wm3 =~ 1*m3
  wy1 =~ 1*y1
  wy2 =~ 1*y2
  wy3 =~ 1*y3

  wx2 ~ ax*wx1 + mx*wm1 + yx*wy1
  wm2 ~ xm*wx1 + am*wm1 + ym*wy1
  wy2 ~ xy*wx1 + my*wm1 + ay*wy1

  wx3 ~ ax*wx2 + mx*wm2 + yx*wy2
  wm3 ~ xm*wx2 + am*wm2 + ym*wy2
  wy3 ~ xy*wx2 + my*wm2 + ay*wy2

  wx1 ~~ wx1
  wx2 ~~ vx*wx2
  wx3 ~~ vx*wx3
  wm1 ~~ wm1
  wm2 ~~ vm*wm2
  wm3 ~~ vm*wm3
  wy1 ~~ wy1
  wy2 ~~ vy*wy2
  wy3 ~~ vy*wy3

  wx1 ~~ wm1
  wx1 ~~ wy1
  wm1 ~~ wy1
  wx2 ~~ cxm2*wm2
  wx2 ~~ cxy2*wy2
  wm2 ~~ cmy2*wy2
  wx3 ~~ cxm3*wm3
  wx3 ~~ cxy3*wy3
  wm3 ~~ cmy3*wy3
'

##########################################################
### Model 4: innovation covariances constrained equal    ###
##########################################################

riclpm_model4 <- '
  RIx =~ 1*x1 + 1*x2 + 1*x3
  RIm =~ 1*m1 + 1*m2 + 1*m3
  RIy =~ 1*y1 + 1*y2 + 1*y3

  x1 ~ mux1*1
  x2 ~ mux2*1
  x3 ~ mux3*1
  m1 ~ mum1*1
  m2 ~ mum2*1
  m3 ~ mum3*1
  y1 ~ muy1*1
  y2 ~ muy2*1
  y3 ~ muy3*1

  RIx ~~ RIx
  RIm ~~ RIm
  RIy ~~ RIy
  RIx ~~ RIm
  RIx ~~ RIy
  RIm ~~ RIy

  wx1 =~ 1*x1
  wx2 =~ 1*x2
  wx3 =~ 1*x3
  wm1 =~ 1*m1
  wm2 =~ 1*m2
  wm3 =~ 1*m3
  wy1 =~ 1*y1
  wy2 =~ 1*y2
  wy3 =~ 1*y3

  wx2 ~ ax*wx1 + mx*wm1 + yx*wy1
  wm2 ~ xm*wx1 + am*wm1 + ym*wy1
  wy2 ~ xy*wx1 + my*wm1 + ay*wy1

  wx3 ~ ax*wx2 + mx*wm2 + yx*wy2
  wm3 ~ xm*wx2 + am*wm2 + ym*wy2
  wy3 ~ xy*wx2 + my*wm2 + ay*wy2

  wx1 ~~ wx1
  wx2 ~~ vx*wx2
  wx3 ~~ vx*wx3
  wm1 ~~ wm1
  wm2 ~~ vm*wm2
  wm3 ~~ vm*wm3
  wy1 ~~ wy1
  wy2 ~~ vy*wy2
  wy3 ~~ vy*wy3

  wx1 ~~ wm1
  wx1 ~~ wy1
  wm1 ~~ wy1
  wx2 ~~ cxm*wm2
  wx2 ~~ cxy*wy2
  wm2 ~~ cmy*wy2
  wx3 ~~ cxm*wm3
  wx3 ~~ cxy*wy3
  wm3 ~~ cmy*wy3
'

model1 <- fit_riclpm(riclpm_model1)
model2 <- fit_riclpm(riclpm_model2)
model3 <- fit_riclpm(riclpm_model3)
model4 <- fit_riclpm(riclpm_model4)

fits <- list(
  "Model 1" = model1,
  "Model 2" = model2,
  "Model 3" = model3,
  "Model 4" = model4
)

riclpm_fit <- fit_table(fits)
print(riclpm_fit)


write.xlsx(riclpm_fit, file = "敏感性分析RI_CLPM_fit_comparison.xlsx", rowNames = FALSE)

cat("\nScaled likelihood-ratio tests:\n")
print(lavTestLRT(model1, model2, model3, model4))

cat("\nModel summaries:\n")
summary(model1, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)
summary(model2, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)
summary(model3, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)
summary(model4, standardized = TRUE, fit.measures = TRUE, rsquare = TRUE)

write.csv(
  parameterEstimates(model1, standardized = TRUE, ci = TRUE),
  "RI_CLPM_pos_neg_swb_model1_parameters.csv",
  row.names = FALSE
)
# write.csv(
#   parameterEstimates(model2, standardized = TRUE, ci = TRUE),
#   "RI_CLPM_model2_parameters.csv",
#   row.names = FALSE
# )
write.xlsx(parameterEstimates(model2, standardized = TRUE, ci = TRUE),
           file = "RI_CLPM_model2_parameters.xlsx", rowNames = FALSE)

model2_parameters <- parameterEstimates(model2, standardized = TRUE,ci = TRUE)
write.xlsx( model2_parameters,file = "RI_CLPM_model2_parameters.xlsx",rowNames = FALSE)

write.csv(
  parameterEstimates(model3, standardized = TRUE, ci = TRUE),
  "RI_CLPM_pos_neg_swb_model3_parameters.csv",
  row.names = FALSE
)
write.csv(
  parameterEstimates(model4, standardized = TRUE, ci = TRUE),
  "RI_CLPM_pos_neg_swb_model4_parameters.csv",
  row.names = FALSE
)


###################################################################################   Model 2
pe_model2 <- parameterEstimates(
  model2,
  standardized = TRUE,
  ci = TRUE
)

cross_lagged_model2 <- pe_model2[
  pe_model2$op == "~" &
    (
      (pe_model2$lhs == "wx2" & pe_model2$rhs %in% c("wm1", "wy1")) |
        (pe_model2$lhs == "wm2" & pe_model2$rhs %in% c("wx1", "wy1")) |
        (pe_model2$lhs == "wy2" & pe_model2$rhs %in% c("wx1", "wm1")) |
        (pe_model2$lhs == "wx3" & pe_model2$rhs %in% c("wm2", "wy2")) |
        (pe_model2$lhs == "wm3" & pe_model2$rhs %in% c("wx2", "wy2")) |
        (pe_model2$lhs == "wy3" & pe_model2$rhs %in% c("wx2", "wm2"))
    ),
]

cross_lagged_model2$path <- paste(
  cross_lagged_model2$rhs,
  "->",
  cross_lagged_model2$lhs
)

cross_lagged_model2$p_fdr_bh <- p.adjust(
  cross_lagged_model2$pvalue,
  method = "BH"
)

cross_lagged_model2$significant_fdr_05 <- cross_lagged_model2$p_fdr_bh < 0.05

cross_lagged_fdr_results <- cross_lagged_model2[
  ,
  c(
    "path",
    "lhs",
    "op",
    "rhs",
    "label",
    "est",
    "se",
    "z",
    "pvalue",
    "p_fdr_bh",
    "significant_fdr_05",
    "std.all",
    "ci.lower",
    "ci.upper"
  )
]

print(cross_lagged_fdr_results)

write.csv(
  cross_lagged_fdr_results,
  "RI_CLPM_model2_cross_lagged_FDR_BH.csv",
  row.names = FALSE
)