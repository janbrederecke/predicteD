bl_eventprob <- function(.cox_model
                         , .timepoint = 10
                         , .newdata
                         , .age_var
){

  # Make necessary objects
  surv_curves <- survfit(.cox_model, newdata = .newdata)
  
  n_newdata <- nrow(.newdata)
  surv_age <- surv_age_plus_timepoint <- vector(length = n_newdata)
  
  # Calculate in newdata
  for (i in seq_len(n_newdata)) {
    surv_age[i] <-
      summary(surv_curves, times = .newdata[[.age_var]][i])$surv[, i]
    
    surv_age_plus_timepoint[i] <-
      summary(surv_curves,
              times = .newdata[[.age_var]][i] + .timepoint,
              extend = TRUE)$surv[, i] 
  }
  
  # Return event probability
  1 - surv_age_plus_timepoint / surv_age
}


