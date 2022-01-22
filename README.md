# Programming-assignment: Illustrating problems in the adjustment of covariates when doing causal inference
Authors: Patricia Fernandez Moreno, Sofía González Matatoros, Víctor Manuel López Molina and Álvaro Pita Ramírez

# 1. Introduction to causal inference
According to the Encyclopedia Britannica (2022), induction is a means of reasoning from particulars to generals and one example of this is causal inference. The objective of this process is to infer the actual effect of a particular element within a more complex system, considering the global outcome under changing conditions (Pearl & Judea, 2009). At this point, it is important to clarify that association does not imply causation, in other words, two elements can be related because they belong to the same interacting network, but its relationship cannot be due to a cause-effect process.

Back to the initial topic, when we want to make an interpretation of the statistical results in this type of scenarios, a relevant aspect is to choose on which covariates we are going to adjust in the analysis. Therefore, the aim of this report is to explain which are the consequences of doing an incorrect analysis, since we are aware that the process of causal inference is highly used in biomedical research.

# 2. Explaining our DAG

# 3. Choosing covariates

## 3.1. Example with multiple regression

Given the diagram explained above, we are going to check whether the effect of X and Y over Z is better predicted by fitting or not by the cofounder W. To this end, we will execute the function dist_multiple() described in the following box

```

dist_multiple <- function(N = 50, b_xw = 3, b_yw = 2, b_zy = 5, b_zx = 4, 
                            sd_x = 1, sd_y = 5, sd_z = 2, B = 2000, w_min = 1, 
                            w_max = 5) {
  
  X_no_W <- rep(NA, B) 
  # we initialize the list with the coeficents of Z when the model is 
  # not adjusted
  X_yes_W <- rep(NA, B) 
  # we initialize the list of p-values the coeficents of Z when the model is 
  # not adjusted
  pv_X_no_W <- rep(NA, B) 
  # we initialize the list with the coeficents of Z when the model is adjusted
  pv_X_yes_W <- rep(NA, B) 
  # we initialize the list of p-values the coeficents of Z when the model is
  # adjusted
  
  for(i in 1:B) { 
    W <- runif(N, w_min, w_max) # The variable W follows an uniform distribution
    X <- b_xw * W + rnorm(N, 0, sd = sd_x) 
    # The variable X follows a normal distribution, influenced by W
    Y <- b_yw * W + rnorm(N, 0, sd = sd_y) 
    # The variable Y follows a normal distribution, influenced by W
    Z <- b_zy * Y + b_zx * X + rnorm(N, 0, sd = sd_z)
    # The variable Z follows a normal distribution, 
    # acting as a collider of X and Y
    
    m1 <- lm(Z ~ (X + Y)) 
    # linear regression of Y depending on X and Y without adjusting
    m2 <- lm(Z ~ (X + Y) + W) 
    # linear regression of Y depending on X and Y adjusting by W
    
    X_no_W[i] <- coefficients(m1)["X"] 
    #estimated coefficient of the effect in model m1
    X_yes_W[i] <- coefficients(m2)["X"] 
    #estimated coefficient of the effect in model m1
    
    pv_X_no_W[i] <- summary(m1)$coefficients["X", "Pr(>|t|)"] 
    #p-value in model m1
    pv_X_yes_W[i] <- summary(m2)$coefficients["X", "Pr(>|t|)"] 
    #p-value in model m2
    
    rm(Z, X, Y, W) # We prepare the variables for the next call of the function
  }
  cat("\n Summary effect without W\n")
  print(summary(X_no_W)) #summary of the coefficents without adjust by W
  cat("\n s.d. estimate = ", sd(X_no_W)) #standard deviation without adjustment
  cat("\n\n Summary effect with W\n")
  print(summary(X_yes_W)) #summary of the coefficents without adjust by W
  cat("\n s.d. estimate = ", sd(X_yes_W), "\n") 
  #standard deviation with adjustment
  
  op <- par(mfrow = c(2, 2))
  hist(X_no_W, main = "Effect without W in the model", xlab = "Estimate")
  abline(v = b_zx, lty = 2)
  # histogram of the unadjusted coefficents
  hist(X_yes_W, main = "Effect with W in the model", xlab = "Estimate") 
  abline(v = b_zx, lty = 2)
  # histogram of the adjusted coefficents
  
  hist(pv_X_no_W, main = "Effect without W in the model", xlab = "p-value") 
  # histogram of the unadjusted p-values.
  hist(pv_X_yes_W, main = "Effect with W in the model", xlab = "p-value") 
  # histogram of the adjusted p-values.
  par(op)
}

```
When we call the function ``` dist_multiple() ```, the following output is:


### Summary X without W

| Min. | 1st Qu. | Median | Mean | 3rd Qu. | Max.
| --- | ---| --- | --- | --- | --- |
| 3.688 | 3.944 | 4.001 | 3.999 | 4.057 | 4.382 |

s.d. estimate =  0.0866533


### Summary X with W

| Min. | 1st Qu. | Median | Mean | 3rd Qu. | Max.
| --- | ---| --- | --- | --- | --- |
| 2.878 | 3.811 | 4.002 | 4.003 | 4.203 | 5.410 | 

s.d. estimate =  0.2954924 

$INSERTAR FOTO


## 3.1. Example with simple regression

# 4. Changing parameters


# 5. Comparisons and conclusions

# Bibliography
Britannica, T. Information Architects of Encyclopaedia. (2022). thought. *Encyclopedia Britannica*. https://www.britannica.com/facts/thought

Pearl, Judea. (2009). Causal inference in statistics: An overview. *Statistics Surveys*, *3*, 96–146. doi:10.1214/09-SS057.
