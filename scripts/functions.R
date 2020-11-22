## Created by Liz Bowman
## June 30, 2017
## various functions


#=========================================================================================
# Linear regression ------------------
#=========================================================================================

# < Run linear regression on all predictors in data frame > ------------------------------
# data: the data frame to analyze
# response: the column name of the response variable
# 
RegressSimple <- function(data, predictor) {
  #predictor.index <- which(colnames(data) == predictor)
  responses <- colnames(data)[7:length(data)]
  model.list <- list()
  for(response in responses){
    simple <- lm(data[,response] ~ data[,predictor])
    #--extract values of interst (p value, df, r-squared)
    model.list[[response]] <- summary(simple) 
    #cat('Response: ', response, '\n')
  }
  class(model.list) <- 'RegressSimple'
  return(model.list)
}

# < output stats for regressions from RegressSimple object > -----------------------------
print.RegressSimple <- function(x,...) {
  # Get vector of all the elements' names
  responses <- names(x)
  # set up a data frame to store the statistics of interest
  model.results <- data.frame(variable = responses,
                              r2 = NA,
                              p = NA)
  # for(predictor in predictors){
  #   predictor.model<- x[[predictor]]
  #   model.results[model.results$variable == predictor, 'r2'] <- 
  #     round(predictor.model$r.squared, 5)
  #   model.results[model.results$variable == predictor, 'p'] <-
  #     round(predictor.model$coefficients[2,4], 5)
  # }
  # return(model.results)
  
  # extract r-squared and p-values
  model.results$r2 <- sapply(x, '[[', 'r.squared')
  model.coeffs <- sapply(x, '[[','coefficients')
  model.results$p <- model.coeffs[8,]
  return(model.results)
  
  # Print values
  # cat('Regression results:', '\n')
  # print(as.matrix(model.results), quote = F)
}

# sapply returns elements as a vector (1-dimensional) or a matrix (2 dimensional)


#========================================================================================#
# Normtester: test for normality of variables----
#========================================================================================#

normtest <- function(data, factor){
  par(mfrow=c(1,2))
  ## Have a look at the densities
  plot(density(data[,factor]))
  ## Plot using a qqplot
  qqnorm(data[,factor]); qqline(data[,factor], col = 2) # normality
  test <- shapiro.test(data[,factor]) # normality
  print(factor)
  print(test)
}