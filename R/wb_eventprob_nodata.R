#' @title wb_eventprob_nodata
#'
#' @description Retrieves the event probabilites without the original data
#'
#' @param .wb_info The output of a former wb_get_info() function.
#' package.
#' @param .timepoint The desired timepoint for the event probabilities.
#' @param .newdata
#' @param .newdata_timestart
#' @export
#'

wb_eventprob_nodata <- function(.wb_info
                              , .timepoint = 10
                              , .newdata
                              , .newdata_timestart
    ){
    # Add error message in case type is != counting
    mssg <- "The model uses wrong type of censoring"
    if (attr(.wb_info[["fit"]]$y, "type") != "counting")
      stop(mssg)

    # Make necessary objects
    n_newdata <- nrow(.newdata)
    linear_predictor_exp <-
      a_newdata <-
      b_newdata <- lp_vector <- vector(length = n_newdata)
    rm(n_newdata)

    # Calculate linear predictor in newdata
    for (i in 1:nrow(.newdata)) {
      stratum <- .newdata[i, .wb_info[["strata_variable"]]]

      ## Add the Wb parameters by stratum and row in newdata
      a_newdata[i] <-
        .wb_info[["wb_param"]][paste0(.wb_info[["strata_variable"]],
                                     "=",
                                     stratum), ][[1]]
      b_newdata[i] <-
        .wb_info[["wb_param"]][paste0(.wb_info[["strata_variable"]],
                                     "=",
                                     stratum), ][[2]]

      ## Calculate linear predictor per person using strata-means
      for (j in seq_along(.wb_info[["model_coeff"]])) {
        coefficient <- .wb_info[["model_coeff"]][j]
        mean <-
          unlist(.wb_info[["means"]][[paste0(.wb_info[["strata_variable"]],
                                            "=",
                                            stratum)]][j])
        value <-
          .newdata[i, colnames(.wb_info[["model_coeff"]])[j]]

        lp_part <- (value - mean) * coefficient

        lp_vector[j] <- lp_part
      }
      linear_predictor_exp[i] <- exp(sum(lp_vector))
    }
    rm(i, j, stratum, coefficient, mean, value, lp_part, lp_vector)

    # Predict the event-probability
    timestart <- .newdata[, .newdata_timestart]
    1 - exp(- (((
      timestart + .timepoint
    ) / b_newdata) ^ a_newdata -
      (timestart / b_newdata) ^ a_newdata) * linear_predictor_exp)
  }