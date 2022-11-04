#' @title wb_get_info
#'
#' @description Generic function to create HTML tables from lists of results
#'
#' @param .fit
#' @param .timestart
#' @param .timestop
#' @param .status
#' @param .weights
#' @param .inits
#' @param .method
#' @param .control
#' @param .proportional
#' @param .params_greater_one
#' @export
#'

wb_get_info <- function(.fit
                      , .timestart
                      , .timestop
                      , .status
                      , .weights = NULL
                      , .inits = log(c(6, 100))
                      , .method = "Nelder-Mead"
                      , .control = list(fnscale = -1
                                        , .parscale = c(0.1, 0.1)
                                        , .maxit = 2000)
                      , .proportional = TRUE
                      , .params_greater_one = FALSE
    ){

    # Intitialize output-list
    output <- list()

    # Retrieve the model coefficients
    output[["model_coeff"]] <- rbind(.fit[["coefficients"]])

    # Add strata of the original data
    output[["strata"]] <- pmisc::get_strata(.fit)

    # Add strata_variable
    strata_vector <- levels(output[["strata"]])[1]
    output[["strata_variable"]] <-
      unlist(stringr::str_split(strata_vector,
                                pattern = "="))[1]
    rm(strata_vector)

    # Add linear predictor (just to ensure the manual calculation is right)
    output[["lp"]] <- stats::predict(.fit, type = "lp")

    # Add manually calculated lp_manual
    output[["lp_manual"]] <- vector(length = nrow(.fit[["model"]]))
    names(output[["lp_manual"]]) <- rownames(.fit[["model"]])

    ## Calculate means per stratum
    n_levels_strata <- nlevels(output[["strata"]])
    output[["means"]] <- list()

    for (i in 1:n_levels_strata) {
      level <- levels(output[["strata"]])[i]
      v <- which(output[["strata"]] == level)
      data_strata <- .fit$model[v, ]
      output[["means"]][[i]] <-
        lapply(data_strata[, colnames(output[["model_coeff"]])], mean)
      names(output[["means"]])[i] <- level
      rm(v, level, data_strata)
    }

    rm(i)

    lp_vector <- vector(length = length(.fit[["coefficients"]]))

    for (i in 1:nrow(.fit[["model"]])) {
      stratum <- .fit[["model"]][i, paste0("strata(",
                                          output[["strata_variable"]][1],
                                          ")")]

      for (j in seq_along(output[["model_coeff"]])) {
        coefficient <- output[["model_coeff"]][j]
        mean <- unlist(output[["means"]][[stratum]][j])
        value <-
          .fit[["model"]][i, colnames(output[["model_coeff"]])[j]]

        lp_part <- (value - mean) * coefficient

        lp_vector[j] <- lp_part
      }

      output[["lp_manual"]][i] <- sum(lp_vector)
    }

    # Add the Weibull parameters per stratum
    if (is.null(output[["strata"]])) {
      output[["strata"]] <- factor(rep(1, nrow(.fit[["model"]])))
    }

    n_levels <- nlevels(output[["strata"]])

    wb_param <- stratum <- a <- b <- vector()

    for (i in 1:n_levels) {
      v <- which(output[["strata"]] == levels(output[["strata"]])[i])
      wf <- pmisc::Wb_fit(
        timestart = .timestart[v],
        timestop = .timestop[v],
        status = .status[v],
        expLP = exp(output[["lp"]][v]),
        weights = .weights[v],
        inits = .inits,
        method = .method,
        control = .control,
        proportional = .proportional,
        params_greater_one = .params_greater_one
      )
      a[v] <- wf$a
      b[v] <- wf$b
      stratum[v] <- levels(output[["strata"]])[i]

      wb_param <- rbind(wb_param, wf)
      rownames(wb_param)[i] <- levels(output[["strata"]])[i]
    }

    output[["wb_param_sample"]] <- cbind(stratum, a, b)

    # Add the general Weibull parameters for simple use
    output[["wb_param"]] <- wb_param

    # Add a complete parameter matrix per stratum
    output[["param_matrix"]] <- cbind(output[["wb_param"]],
                                      output[["model_coeff"]])

    # Add the model to the output
    output[["fit"]] <- .fit

    # Remove original data from model
    output[["fit"]][["model"]] <- NULL

    # Return output
    output
  }