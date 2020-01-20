cma <- function(edata, outcome, treatment, mediator, result_path_file) {
  
  # An R function that utilizes GMM to estimate causal mediation
  # The inputs are:
  #   edata: experimental unit-level data
  #   outcome: outcome metric of experimental unit
  #   treatment: treatment indicator of experimental unit
  #   mediator: mediator metric of experimental unit
  #   result_path_file: e.g., "~/results/gbdt_desktop_conversion_cma.csv"
  # The output is a CSV file of estimation results

  library(data.table)
  library(lmtest)
  library(sandwich)
  library(gmm)
  
  edata <- edata[, c(outcome, treatment, mediator), with=FALSE]
  
  setnames(edata, c("outcome", "treatment", "mediator"))
  #############################################################################################################################
  ## Descriptive Stats
  des <- edata[, list(mean.outcome = mean(outcome)), by=treatment]
  mean.control = des[treatment==0, mean.outcome]
  cat(paste0(outcome, " Mean in Control Group is ", mean.control), "\n")
  rm(des)
  #############################################################################################################################
  gmm_mediation <- function(delta, d) {
    delta.m0 <- delta[1]
    delta.m1 <- delta[2]
    delta.y0 <- delta[3]
    delta.y1 <- delta[4]
    delta.y2 <- delta[5]
    delta.y3 <- delta[6]
    
    med.moment1 <- d[, mediator] - (delta.m0 * d[, constant]) - (delta.m1 * d[, treatment])
    med.moment2 <- d[, treatment] * (d[, mediator] - (delta.m0 * d[, constant]) - (delta.m1 * d[, treatment]))
    out.moment1 <- d[, outcome] - (delta.y0 * d[, constant]) - (delta.y1 * d[, treatment]) - (delta.y2 * d[, mediator]) - (delta.y3 * d[, mediator] *  d[, treatment])
    out.moment2 <- d[, treatment] * (d[, outcome] - (delta.y0 * d[, constant]) - (delta.y1 * d[, treatment]) - (delta.y2 * d[, mediator]) - (delta.y3 * d[, mediator] *  d[, treatment]))
    out.moment3 <- d[, mediator] * (d[, outcome] - (delta.y0 * d[, constant]) - (delta.y1 * d[, treatment]) - (delta.y2 * d[, mediator]) - (delta.y3 * d[, mediator] *  d[, treatment]))
    out.moment4 <- d[, mediator] * d[, treatment] * (d[, outcome] - (delta.y0 * d[, constant]) - (delta.y1 * d[, treatment]) - (delta.y2 * d[, mediator]) - (delta.y3 * d[, mediator] *  d[, treatment]))
    
    g <- cbind(med.moment1, med.moment2, out.moment1, out.moment2, out.moment3, out.moment4)
    return(g)
  }
  
  #############################################################################################################################
  results_mediator <- lm(mediator ~ treatment, data = edata)
  
  delta.m0.e = unname(results_mediator$coefficients['(Intercept)'])
  delta.m1.e = unname(results_mediator$coefficients['treatment'])
  
  #############################################################################################################################
  results <- lm(outcome ~ treatment + mediator + treatment:mediator, data = edata)
  
  delta.y0.e = unname(results$coefficients['(Intercept)'])
  delta.y1.e = unname(results$coefficients['treatment'])
  delta.y2.e = unname(results$coefficients['mediator'])
  delta.y3.e = unname(results$coefficients['treatment:mediator'])
  #############################################################################################################################
  
  edata <- edata[, constant:=1.000]
  
  start_time <- Sys.time()
  results <- gmm(gmm_mediation, 
                 edata, 
                 c(delta.m0 = delta.m0.e, 
                   delta.m1 = delta.m1.e, 
                   delta.y0 = delta.y0.e, 
                   delta.y1 = delta.y1.e, 
                   delta.y2 = delta.y2.e, 
                   delta.y3 = delta.y3.e
                 ), 
                 traceIter = TRUE,
                 wmatrix = "optimal",
                 type = "twoStep", 
                 vcov = "HAC")
  end_time <- Sys.time()
  
  rm(edata)
  
  duration = end_time - start_time
  convergence = results$algoInfo$convergence
  
  print(summary(results))
  #############################################################################################################################
  delta.m0.e = unname(results$coefficients['delta.m0'])
  delta.m1.e = unname(results$coefficients['delta.m1'])
  
  delta.y0.e = unname(results$coefficients['delta.y0'])
  delta.y1.e = unname(results$coefficients['delta.y1'])
  delta.y2.e = unname(results$coefficients['delta.y2'])
  delta.y3.e = unname(results$coefficients['delta.y3'])
  
  var.m0 = results$vcov['delta.m0', 'delta.m0']
  var.m1 = results$vcov['delta.m1', 'delta.m1']
  var.y0 = results$vcov['delta.y0', 'delta.y0']
  var.y1 = results$vcov['delta.y1', 'delta.y1']
  var.y2 = results$vcov['delta.y2', 'delta.y2']
  var.y3 = results$vcov['delta.y3', 'delta.y3']
  
  covar.m0.m1 = results$vcov['delta.m0', 'delta.m1']
  covar.m0.y0 = results$vcov['delta.m0', 'delta.y0']
  covar.m0.y1 = results$vcov['delta.m0', 'delta.y1']
  covar.m0.y2 = results$vcov['delta.m0', 'delta.y2']
  covar.m0.y3 = results$vcov['delta.m0', 'delta.y3']
  
  covar.m1.y0 = results$vcov['delta.m1', 'delta.y0']
  covar.m1.y1 = results$vcov['delta.m1', 'delta.y1']
  covar.m1.y2 = results$vcov['delta.m1', 'delta.y2']
  covar.m1.y3 = results$vcov['delta.m1', 'delta.y3']
  
  covar.y0.y1 = results$vcov['delta.y0', 'delta.y1']
  covar.y0.y2 = results$vcov['delta.y0', 'delta.y2']
  covar.y0.y3 = results$vcov['delta.y0', 'delta.y3']
  
  covar.y1.y2 = results$vcov['delta.y1', 'delta.y2']
  covar.y1.y3 = results$vcov['delta.y1', 'delta.y3']
  
  covar.y2.y3 = results$vcov['delta.y2', 'delta.y3']
  
  #############################################################################################################################
  
  mediation_result <- data.table(matrix(double(), 5, 3))
  colnames(mediation_result) <- c("Effect",
                                  "Estimate",
                                  "SE")
  
  mediation_result <- mediation_result[, Effect:=as.character(Effect)]
  set(mediation_result, i = 1L, j = "Effect", value =  "ADE_control")
  set(mediation_result, i = 2L, j = "Effect", value =  "ADE_treated")
  set(mediation_result, i = 3L, j = "Effect", value =  "AME_control")
  set(mediation_result, i = 4L, j = "Effect", value =  "AME_treated")
  set(mediation_result, i = 5L, j = "Effect", value =  "Total Effect")
  
  set(mediation_result, i = 1L, j = "Estimate", value = delta.y1.e + (delta.y3.e * (delta.m0.e + (delta.m1.e * 0))))
  set(mediation_result, i = 2L, j = "Estimate", value = delta.y1.e + (delta.y3.e * (delta.m0.e + (delta.m1.e * 1))))
  set(mediation_result, i = 3L, j = "Estimate", value = delta.m1.e * (delta.y2.e + (delta.y3.e * 0)))
  set(mediation_result, i = 4L, j = "Estimate", value = delta.m1.e * (delta.y2.e + (delta.y3.e * 1)))
  
  #############################################################################################################################
  
  var.ame <- function(t) {
    return((((delta.y2.e + (delta.y3.e * t))^2) * var.m1)
           + ((delta.m1.e^2) * var.y2)
           + (((delta.m1.e * t)^2) * var.y3)
           + (2 * (delta.y2.e + (delta.y3.e * t)) * delta.m1.e * covar.m1.y2)
           + (2 * (delta.y2.e + (delta.y3.e * t)) * delta.m1.e * t * covar.m1.y3)
           + (2 * (delta.m1.e^2) * t * covar.y2.y3)
    )
  }
  
  
  var.ade <- function(t) {
    return(var.y1
           + (((delta.m0.e + (delta.m1.e * t))^2) * var.y3)
           + ((delta.y3.e^2) * var.m0)
           + (((delta.y3.e * t)^2) * var.m1)
           + (2 * delta.y3.e * covar.m0.y1)
           + (2 * (delta.m0.e + (delta.m1.e * t)) * covar.y1.y3)
           + (2 * delta.y3.e * t * covar.m1.y1)
           + (2 * (delta.m0.e + (delta.m1.e * t)) * delta.y3.e * covar.m0.y3)
           + (2 * (delta.m0.e + (delta.m1.e * t)) * delta.y3.e * t * covar.m1.y3)
           + (2 * (delta.y3.e^2) * t * covar.m0.m1)
    )
  }
  
  set(mediation_result, i = 1L, j = "SE", value = sqrt(var.ade(0)))
  set(mediation_result, i = 2L, j = "SE", value = sqrt(var.ade(1)))
  set(mediation_result, i = 3L, j = "SE", value = sqrt(var.ame(0)))
  set(mediation_result, i = 4L, j = "SE", value = sqrt(var.ame(1)))
    
  mediation_result[, "Z_Score":= Estimate/SE]  
  mediation_result[, "P_Value":= 2*pnorm(-abs(Z_Score))]
  mediation_result[, "% Change":= Estimate/mean.control]
  
  #############################################################################################################################
  
  lm <- lm(outcome ~ treatment, data = results$dat)
  ttest <- coeftest(lm, vcov = vcovHAC(lm))
  
  set(mediation_result, i = 5L, j = "Estimate", value = ttest["treatment", "Estimate"])
  set(mediation_result, i = 5L, j = "SE", value = ttest["treatment", "Std. Error"])
  set(mediation_result, i = 5L, j = "Z_Score", value = ttest["treatment", "t value"])
  set(mediation_result, i = 5L, j = "P_Value", value = ttest["treatment", "Pr(>|t|)"])
  set(mediation_result, i = 5L, j = "% Change", value = ttest["treatment", "Estimate"]/mean.control)
  
  #############################################################################################################################
  cat(paste0("The outcome metric is ", outcome), "\n")
  cat(paste0("The mediating metric is ", mediator), "\n")  
  cat("Mediation Results:", "\n")
  print(mediation_result)
  
  cat(paste0("Convergence Code is ", convergence), "\n")
  if (convergence == 0) {
    cat("The GMM algorithm converged", "\n")
  } else if (convergence == 1) {
    cat("The GMM algorithm did not converge", "\n")
  }

  cat(paste0("Duration is ", duration), "\n")
  
  for (j in 2:dim(mediation_result)[2]) set(mediation_result, j=j, value=round(mediation_result[[j]], 6))

  fwrite(mediation_result, result_path_file)
  
  return(mediation_result)
}
