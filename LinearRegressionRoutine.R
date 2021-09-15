# To quickly re-acquaint myself with some basics in R, I'll write a function that does linear regression

reg <- function(dv,iv,intercept=TRUE) {

  # Store parameters N and K from the dimensions of the iv object
  N <- dim(iv)[1]
  K <- dim(iv)[2]
  
  # If Intercept, join an array of 1's to IV matrix when intercept=TRUE  
  if (intercept) {
    constant <- 1
    iv <- cbind(iv, constant)
  }

  # Solve for the beta-coefficients,
  betas <- solve(t(iv) %*% iv) %*% (t(iv) %*% dv)

  # Calculate fitted values and residuals
  fitted <- iv %*% betas
  resid <- dv - fitted

  # Use residuals to calculate R-squared and other fit-statistics
  sq_resid <- resid^2
  SSR <- sum(sq_resid)
  SSM <- sum((fitted-mean(dv))^2)
  TSS <- sum((dv-mean(dv))^2)
  
  Rsquared <- 1-(SSR/TSS)
  Adj_Rsquared <- 1 - (((1-Rsquared)*(N-1))/(N-K-1))
  RMSE <- (SSR/(N-K-1))^0.5
  
  F_Statistic <- (SSM/(K))/(SSR/(N-K-1))
  Prob_F <- pf(F_Statistic,K,(N-K-1),lower.tail = FALSE)
  
  # Estimate standard errors, T-stats, and p-values
  sigma_sq <- SSR/(N-K-1)
  varcov <- sigma_sq * (solve(t(iv) %*% iv))
  std_err <- diag(varcov)^0.5
  t_stats <- betas / std_err
  pvals <- 2*pt(abs(t_stats),(N-K-1),lower.tail=FALSE)
  
  
  # Return a list of three table-objects: betas, fit-statics, and fitted-values/residuals
  summary <- cbind(betas,std_err,t_stats,pvals)
  colnames(summary) <- cbind("Betas","S.E.","T-Stats","P-Value")
  
  fit_statistics <- rbind(N,Rsquared,Adj_Rsquared,RMSE,F_Statistic,Prob_F)
  colnames(fit_statistics) <- c("Model Stats")
  
  data <- data.frame(dv,fitted,resid)
  colnames(data) <- c("Actual","Fitted","Residual")
  
  results <- list(summary = summary,fit_statistics = fit_statistics,values = data)
  return(results)
}

# Do a test of my function with some sample data
library(MASS)
df <- Cars93

# First, let's try it with the lm function in R
lm_result <- lm(Price~MPG.city+MPG.highway, data=df)
summary(lm_result)

# Second, let's see how my function performs
my_result <- reg(df$Price,data.matrix(df[c("MPG.city","MPG.highway")]))
my_result$fit_statistics
my_result$summary

