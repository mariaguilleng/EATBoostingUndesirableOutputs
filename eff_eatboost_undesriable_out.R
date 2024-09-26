#' @title Efficiency score for undesriable outputs
#'
#' @description This function calculates the efficiency score in the presence
#' of undesirable outputs assuming in a strong disposability context.
#'
#' @name DDF_T1
#'
#' @param data \code{data.frame} containing the DMUs' data
#' @param x Column of good input indexes in \code{data}.
#' @param y_desirable Column of desirable outputs indexes in \code{data}.
#' @param y_undesirable Column of un desriable outputs indexes in \code{data}.
#' @param model A \code{EATBoost} model from package \code{boostingDEA}.
#'
#' @return The efficiency score for each DMU


eff_eatboost_undesriable_out <- function(data, x, y_desirable, y_undesirable = NULL, model = NULL) {
  
  if (class(model) != "EATBoost") {
    stop("Model must be EATBoost")
  }
  
  # Handling the eatboosting scenario
  xOriginal <- model[["data"]][["x"]]
  yOriginal <- model[["data"]][["y"]]
  dataOriginal <- model[["data"]][["df"]]
  pred <- predict(model, dataOriginal, xOriginal)
  dataOriginal <- cbind(dataOriginal[, xOriginal], pred)
  y_desirable_pred_k <- as.matrix(dataOriginal[, y_desirable])
  
  j <- nrow(dataOriginal)
  xOriginal_k <- as.matrix(dataOriginal[, xOriginal])
  yOriginal_k <- as.matrix(dataOriginal[, yOriginal])
  
  # Handling current scenario
  numDMU <- nrow(data)
  x_k <- as.matrix(data[, x])
  y_desirable_k <- as.matrix(data[, y_desirable])
  scores <- matrix(nrow = numDMU, ncol = 1)
  nX <- length(x)
  nY <- if (is.null(y_undesirable)) length(y_desirable) else length(y_desirable) + length(y_undesirable)
  
  if (!is.null(y_undesirable)) {
    y_undesirable_k <- as.matrix(data[, y_undesirable])
  }
  
  # Get scores
  for (d in 1:numDMU) {
    
    # Objetive function
    objVal <- matrix(ncol = j + 3, nrow = 1)
    objVal[1] <- 1
    objVal[2] <- 0.000001
    objVal[3] <- 0.000001
    
    # Structure for lpSolve
    lps <- make.lp(nrow = nX + nY, ncol = j + 3)
    lp.control(lps, sense = "max")
    set.objfn(lps, objVal)
    
    # Constrain inputs
    for (xi in x) {
      add.constraint(lps, xt = c(0, 1, 0, x_k[, xi]), "=", rhs = x_k[d, xi])
    }
    
    # Constraint undesirable outputs 
    if (!is.null(y_undesirable)) {
      for (yi in 1:length(y_undesirable)) {
        add.constraint(lps, xt = c(y_undesirable_k[d, yi], 0, 0, y_undesirable_k[, yi]), "=", rhs = 2 * y_undesirable_k[d, yi])
      }
    }
    
    # Constrain desirable outputs
    for (yi in 1:length(y_desirable)) {
      add.constraint(lps, xt = c(-y_desirable_k[d, yi], 0, -1, y_desirable_pred_k[, yi]), "=", rhs = 0)
    }
    
    # Constrain sum(lambda) = 1
    add.constraint(lps, xt = c(0, 0, 0, rep(1, j)), type = "=", rhs = 1)
    
    # Solution
    errorCode <- solve(lps)
    scores[d, ] <- round(get.objective(lps), 4)
    write.lp(lps, "new.lp")
  }
  
  return(scores)
}
