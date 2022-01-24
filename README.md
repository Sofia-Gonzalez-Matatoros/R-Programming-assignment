# Programming-assignment: Illustrating problems in the adjustment of covariates when doing causal inference
Authors: Patricia Fernandez Moreno, Sofía González Matatoros, Víctor Manuel López Molina and Álvaro Pita Ramírez

# 1. Introduction 

## 1.1. What is causal inference?

According to the Encyclopedia Britannica (2022), induction is a means of reasoning from particulars to generals and one example of this is causal inference. The objective of this process is to infer the actual effect of a particular element within a more complex system, considering the global outcome under changing conditions (Pearl & Judea, 2009). At this point, it is important to clarify that association does not imply causation, in other words, two elements can be related because they belong to the same interacting network, but its relationship cannot be due to a cause-effect process.

## 1.2. Covariates: Fork, collider and backdoor path

Before starting with the objectives, it is essential to understand the following concepts clearly, since they will appear throughout the report.

- Fork: diagram in which two arrows emanate from a variable to two independent nodes, giving rise to a dichotomy. (Pearl et al., 2016)
- Collider: structure in which one variable receives edges from two other nodes. (Pearl et al., 2016)
- Backdoor path: given a graph, a backdoor path from one point A to another B is a non-causal path that begins with a parent of A and ends at B with the collaboration of other points that act as mediators. If we want to identify the direct effect (causal path) of A over B when existing backdoor paths, the backdoor criterion states that the rest of measured covariates must be sufficient to block all backdoor paths from A to B. (Blackwell, 2013)

## 1.3. Objective

Back to the initial topic, when we want to make an interpretation of the statistical results in this type of scenarios, a relevant aspect is to choose on which covariates we are going to adjust in the analysis. Therefore, the aim of this report is to explain which are the consequences of doing an incorrect analysis, since we are aware that the process of causal inference is highly used in biomedical research.


# 2. Explaining our DAG

To illustrate how to perform the correct analysis, we are going to examinate the following scenario described on this code

```
library(dagitty)
if(!suppressWarnings(require("rethinking", quietly = TRUE))) {
  drawdag <- plot
} 

dag1 <- dagitty("dag {
W -> X
W -> Y
X -> Z
Y -> Z
e_x -> X
e_y -> Y
e_zx -> Z
e_zy -> Z
}")

coordinates(dag1) <- list(x = c(W = 2.5, X = 2, Y = 3, Z = 2.5, e_x = 1.5, 
                                e_y = 3.5, e_zx =1.5, e_zy = 3.5),
                          y = c(W = -0.5, X = 0, Y = 0, Z = 0.5, e_x = -.25, 
                                e_y = -0.25, e_zx = .5, e_zy = .5))
## plot(dag1)
drawdag(dag1)
drawdag(dag1, xlim = c(1, 4), ylim = c(-1, 1))

```


![Rplot03](https://user-images.githubusercontent.com/97886286/150672937-08e63785-6c48-4361-91a0-9c9ebad09cc4.jpeg)

Considering this DAG, we can make these observations
- X and Y are associated throw W, which is a fork, so conditioning on W would render them independent.
- Z is a collider, depending on X and Y.

# 3. Choosing covariates

## 3.1. Example with multiple regression

Given the diagram explained above, we are going to check whether the effect of X and Y over Z is better predicted by fitting or not by the cofounder W. To this end, we will execute the function dist_multiple() described in the following box

```

dist_multiple <- function(N = 50, b_xw = 3, b_yw = 2, b_zy = 5, b_zx = 4, 
                          sd_x = 1, sd_y = 5, sd_z = 2, B = 2000, w_min = 1, 
                          w_max = 5) {
  
  X_no_W <- rep(NA, B) 
  # we initialize the list with the coefficents of X when the model is 
  # not adjusted
  X_yes_W <- rep(NA, B)
  # we initialize the list with the coefficents of X when the model is adjusted
  pv_X_no_W <- rep(NA, B) 
  # we initialize the list of p-values for the estimations of the coefficents 
  #of X when the model is not adjusted by W
  pv_X_yes_W <- rep(NA, B) 
  # we initialize the list of p-values for the estimations of the coefficents 
  #of X when the model is adjusted by W
  
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
    # linear regression of Z depending on X and Y without adjusting
    m2 <- lm(Z ~ (X + Y) + W) 
    # linear regression of Z depending on X and Y adjusting by W
    
    X_no_W[i] <- coefficients(m1)["X"] 
    #estimated coefficient of the effect in model m1
    X_yes_W[i] <- coefficients(m2)["X"] 
    #estimated coefficient of the effect in model m2
    
    pv_X_no_W[i] <- summary(m1)$coefficients["X", "Pr(>|t|)"] 
    #p-value in model m1
    pv_X_yes_W[i] <- summary(m2)$coefficients["X", "Pr(>|t|)"] 
    #p-value in model m2
    
    rm(Z, X, Y, W) # We prepare the variables for the next call of the function
  }
  cat("\n Summary effect without W\n")
  print(summary(X_no_W)) #summary of the coefficents without adjusting by W
  cat("\n s.d. estimate = ", sd(X_no_W)) #standard deviation without adjustment
  cat("\n\n Summary effect with W\n")
  print(summary(X_yes_W)) #summary of the coefficents adjusting by W
  cat("\n s.d. estimate = ", sd(X_yes_W), "\n") 
  #standard deviation with the adjustment
  
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


### Summary Effect without W

| Min. | 1st Qu. | Median | Mean | 3rd Qu. | Max.
| --- | ---| --- | --- | --- | --- |
| 3.688 | 3.944 | 4.001 | 3.999 | 4.057 | 4.382 |

s.d. estimate =  0.0866533


### Summary Effect with W

| Min. | 1st Qu. | Median | Mean | 3rd Qu. | Max.
| --- | ---| --- | --- | --- | --- |
| 2.878 | 3.811 | 4.002 | 4.003 | 4.203 | 5.410 | 

s.d. estimate =  0.2954924 

![Plots_multiple](https://user-images.githubusercontent.com/97886286/150673370-9aeefe72-6884-4705-bdd1-2cb0f71c010d.jpeg)


It is observed that when adjusting for W the standard deviation increases 0.11 units, while the mean hardly differs. 
The R documentation indicates that adjusting by a confounder improves the estimates, so we were struck by the fact that in this case adjusting by W increased the deviation. Our hypothesis is that acting W as a fork is part of the information flow from X and Y to Z, since it serves as a mediator between X and Y (Foraita et al., 2014), therefore, adjusting for W would imply eliminating a real association between X and Y, which would affect the information flow received by Z. 


## 3.2. Example with simple regression

In this section we will try to find out the individual effect exerted by X on Z. To this end, we will apply the function ``` dist_simple ```, to see which scenario provides the correct analysis. Here we have four options for adjusting: none, W, Y, or W and Y.

```
dist_simple <- function(N = 50, b_xw = 3, b_yw = 2, b_zy = 5, b_zx = 4, 
                        sd_x = 1, sd_y = 5, sd_z = 2, B = 2000, w_min = 1, 
                        w_max = 5) {
  X_no_W <- rep(NA, B) 
  X_yes_W <- rep(NA, B)
  X_yes_Y <- rep(NA, B)
  X_yes_WY <- rep(NA, B)
  pv_X_no_W <- rep(NA, B) 
  pv_X_yes_W <- rep(NA, B)
  pv_X_yes_Y <- rep(NA, B)
  pv_X_yes_WY <- rep(NA, B)

  
  for(i in 1:B) { 
    W <- runif(N, w_min, w_max)
    X <- b_xw * W + rnorm(N, 0, sd = sd_x)
    Y <- b_yw * W + rnorm(N, 0, sd = sd_y)
    Z <- b_zy * Y + b_zx * X + rnorm(N, 0, sd = sd_z)
    
    m1 <- lm(Z ~ X) # linear regression of Z depending on X (unadjusted)
    m2 <- lm(Z ~ X + W) # adjusted by W
    m3 <- lm(Z ~ X + W + Y) # adjusted by W and Y
    m4 <- lm(Z ~ X + Y) # adjusted by Y
    
    X_no_W[i] <- coefficients(m1)["X"] 
    X_yes_W[i] <- coefficients(m2)["X"] 
    X_yes_WY[i] <- coefficients(m3)["X"]
    X_yes_Y[i] <- coefficients(m4)["X"]
    
    pv_X_no_W[i] <- summary(m1)$coefficients["X", "Pr(>|t|)"] 
    pv_X_yes_W[i] <- summary(m2)$coefficients["X", "Pr(>|t|)"] 
    pv_X_yes_WY[i] <- summary(m3)$coefficients["X", "Pr(>|t|)"] 
    pv_X_yes_Y[i] <- summary(m4)$coefficients["X", "Pr(>|t|)"] 
    
    rm(Z, X, Y, W)
  }
  cat("\n Summary X without adjusting\n")
  print(summary(X_no_W))
  cat("\n s.d. estimate = ", sd(X_no_W))
  cat("\n\n Summary X with W\n")
  print(summary(X_yes_W))
  cat("\n s.d. estimate = ", sd(X_yes_W), "\n")
  cat("\n\n Summary X with W and Y\n")
  print(summary(X_yes_WY))
  cat("\n s.d. estimate = ", sd(X_yes_WY), "\n")
  cat("\n\n Summary X with Y\n")
  print(summary(X_yes_Y))
  cat("\n s.d. estimate = ", sd(X_yes_Y), "\n")
  
  op <- par(mfrow = c(2, 4))
  hist(X_no_W, main = "X without W and Y in the model", xlab = "Estimate")
  abline(v = b_zx, lty = 2)
  hist(X_yes_W, main = "X with W in the model", xlab = "Estimate") 
  abline(v = b_zx, lty = 2)
  hist(X_yes_WY, main = "X with W and Y in the model", xlab = "Estimate") 
  abline(v = b_zx, lty = 2)
  hist(X_yes_Y, main = "X with Y in the model", xlab = "Estimate") 
  abline(v = b_zx, lty = 2)
  
  hist(pv_X_no_W, main = "X without W and Y in the model", xlab = "p-value")
  hist(pv_X_yes_W, main = "X with W in the model", xlab = "p-value")
  hist(pv_X_yes_WY, main = "X with W and Y in the model", xlab = "p-value")
  hist(pv_X_yes_Y, main = "X with Y in the model", xlab = "p-value")
  par(op)
}


```
When we call the function ``` dist_simple() ```, the following output is:

### Summary X without adjusting

| Min. | 1st Qu. | Median | Mean | 3rd Qu. | Max.
| --- | ---| --- | --- | --- | --- |
| 3.820 | 6.404 | 7.098 | 7.082 | 7.775 | 10.655 

s.d. estimate =  1.022683

### Summary X with W

| Min. | 1st Qu. | Median | Mean | 3rd Qu. | Max.
| --- | ---| --- | --- | --- | --- |
|-10.210 | 1.508 | 3.951 | 3.970 | 6.333 | 17.020 

 s.d. estimate =  3.604197 

### Summary X with W and Y

| Min. | 1st Qu. | Median | Mean | 3rd Qu. | Max.
| --- | ---| --- | --- | --- | --- |
| 3.115 | 3.798 | 3.996 | 3.996 | 4.198 | 5.184 

 s.d. estimate =  0.2931457 


### Summary X with Y

| Min. | 1st Qu. | Median | Mean | 3rd Qu. | Max.
| --- | ---| --- | --- | --- | --- |
| 3.679 | 3.941 | 3.999 | 4.000 | 4.057 | 4.365 

 s.d. estimate =  0.08639647


![Plots_simple](https://user-images.githubusercontent.com/97886286/150672678-7e51b079-7291-4dec-98ec-f6015451ded3.jpeg)


Looking at the results, the mean is found to be closer to reality (mean = 4 and s.d. = 0.086) when adjusting by Y. However, comparing the unadjusted model and the models adjusted for at least one covariate, the latter provide more accurate estimates because it blocks the backdoor path (W -> Y -> X). This example is consistent with what has been published in other causal inference articles, in fact, a similar example appears in the documentation of R stored in CRAN (https://cran.r-project.org/web/packages/ggdag/vignettes/intro-to-dags.html) 

# 4. Changing parameters

Considering the examples illustrated above, we are going to check what happens if we decrease the coefficient b_zw and its standard deviation calling the functions like this ``` dist_simple(sd_x = 5, b_zx = -4) ``` and ``` dist_multiple(sd_x = 5, b_zx = -4) ```.

## 4.1. Multiple regression

![Plots_multiple_parametros_cambiados](https://user-images.githubusercontent.com/97886286/150673636-9fd4d2d3-acd2-438e-8b24-03222617dcf8.jpeg)

### Summary effect without W

| Min. | 1st Qu. | Median | Mean | 3rd Qu. | Max.
| --- | ---| --- | --- | --- | --- | 
|-4.162 | -4.034 | -3.999 | -4.000 | -3.966 | -3.835 

 s.d. estimate =  0.05064404

### Summary effect with W

| Min. | 1st Qu. | Median | Mean | 3rd Qu. | Max.
| --- | ---| --- | --- | --- | --- |
| -4.191 | -4.040 | -4.000 | -4.000 | -3.960 | -3.807 

 s.d. estimate =  0.06008477


## 4.2. Simple regression

![Plots_simple_parametros_cambiados](https://user-images.githubusercontent.com/97886286/150673717-eca41ea7-01c8-41d6-92a2-1737e74ddc4c.jpeg)

### Summary X without adjusting

| Min. | 1st Qu. | Median | Mean | 3rd Qu. | Max.
| --- | ---| --- | --- | --- | --- |
| -5.0230 | -3.3424 | -2.9283 | -2.9141 | -2.4988 | -0.2526 

 s.d. estimate =  0.6309767

### Summary X with W

| Min. | 1st Qu. | Median | Mean | 3rd Qu. | Max.
| --- | ---| --- | --- | --- | --- |
| -6.642 | -4.498 | -4.009 | -4.005 | -3.501 | -1.280 

 s.d. estimate =  0.7356117 

### Summary X with W and Y

| Min. | 1st Qu. | Median | Mean | 3rd Qu. | Max.
| --- | ---| --- | --- | --- | --- |
| -4.189 | -4.042 | -4.000 | -4.000 | -3.961 | -3.775 

 s.d. estimate =  0.06130837 

### Summary X with Y

| Min. | 1st Qu. | Median | Mean | 3rd Qu. | Max.
| --- | ---| --- | --- | --- | --- |
| -4.203 | -4.034 | -4.000 | -4.000 | -3.968 | -3.855


 s.d. estimate =  0.05004781
 
## 4.3. Comparison between simple and multiple linear regression

By decreasing the value of the parameter b_zx and increasing the value of sd_x, it is found that the estimates of the deviation of all fitted models improve considerably.

To find out what caused this phenomenon, we reviewed the properties of normal distributions. X and Y are variables that follow a normal distribution:

![image](https://user-images.githubusercontent.com/97886286/150767156-5eb25d8a-56d2-4cb5-aab7-d02efd53886c.png)

If we keep the parameters of Y fixed (mean and deviation) and the deviation of X is increased, we would find that the radicand increases and, therefore, the deviation of the distribution generated from the sum of X and Y. 

Our hypothesis is that the higher the deviation sd_x, the greater the dispersion of the data of the X distribution. Then, considering the causal relationship relating X to Z, the distribution of the data of Z increased the greater the dispersion of X. However, the lower the coefficient b_zx within the natural numbers, the less the influence of X on Y and, therefore, the less the distribution of Z deviates even though the deviation of X is high. Similarly, within negative numbers, the higher the value of b_zx, the lower the influence of X on Y, so the less the distribution of Z would deviate. 

# 5. Conclusions

# 6. Difficulties

In this work, we also wanted to test the adjustment of covariates when doing a logistic regression, to wit, a statistical model in which variables are categorical, such as pass/fail, win/lose, healthy/sick, etc. In this case, we have based our work in the article form Sjölander (2018).

In this article, the author uses logistic regression models to estimate causal effect measures, by the use of the R-package AF. For this purpose, the author, as well as us, uses a publicly available dataset (clslowbwt), which includes information on 487 births among 188 women. Among this information, we can find the variables lbw (a binary indicator of whether the newborn child has low birthweight, defined as a birthweight smaller or equal to 2500 g), smoker (a binary indicator of whether the mother smoked during pregnancy), race (race of the mother, coded as 1. White, 2. Black or 3. Other), age (age of the mother), and id (a unique identification number for each mother). The full information contained in this dataset can be seen in the following link:
https://rdrr.io/cran/AF/man/clslowbwt.html#heading-1 

In our case, we wanted to estimate the proportion of low birth weights that would be prevented if nobody would smoke during pregnancy. We controlled for the mother's race and age in the analysis.

The first step was to install the R-package AF and load it, as well as the clslowbwt dataset we will use for the generation of these models:

```

install.packages("AF")
library(AF)
data(clslowbwt)

```

Once the dataset and the package were loaded, the next step was to fit a logistic regression model that relates the outcome (low birthweight) to the exposure (smoking). This is done by:

```

Simple logistic regression model
model<-glm(formula=lbw~smoker,family="binomial",data=clslowbwt)
summary(model)

```

The result of this is the following:

> summary(model)

Call:
glm(formula = lbw ~ smoker, family = "binomial", data = clslowbwt)

Deviance Residuals:

| Min | 1Q | Median | 3Q | Max
| --- | ---| --- | --- | --- |
-1.0361 | -0.7404 | -0.7404 | 1.3256 | 1.6901  

Coefficients:

| Coefficients |Estimate | Std. Error| z value | Pr(>z)
| --- | ---| --- | --- | --- |
| (Intercept) | -1.1542 | 0.1371 | -8.420 | < 2e-16 *** |
smoker | 0.8124 | 0.1998 | 4.067 | 4.77e-05 ***

Signif. codes:  
0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for binomial family taken to be 1)

Null deviance: 603.05  on 486  degrees of freedom
Residual deviance: 586.36  on 485  degrees of freedom
AIC: 590.36

Number of Fisher Scoring iterations: 4

As we can see, there is strong evidence of the effect of the variable smoker on the outcome (low birthweight). However, we should find out if other variables are also associated with low birthweight. For instance, we will control by race and age:

```

Multiple logistic regression (adjusting by race and age)
model1 <- glm(formula=lbw~smoker+age+race,family="binomial",data=clslowbwt)
summary(model1)

```

The result of this is the following:

summary(model1)

Call:
glm(formula = lbw ~ smoker + age + race, family = "binomial", 
    data = clslowbwt)

Deviance Residuals: 

| Min | 1Q | Median | 3Q | Max
| --- | ---| --- | --- | --- |
| -1.2326 | -0.8936 | -0.6491 | 1.2249 | 1.9808  

Coefficients:

| Coefficients |Estimate | Std. Error| z value | Pr(>z) |
| --- | ---| --- | --- | --- |
| (Intercept) | -1.35946 | 0.53281 | -2.551 | 0.01073 * |
| smoker | 0.60080 | 0.21693 | 2.770 | 0.00561 **
age | 0.02399 | 0.01785 | 1.344 | 0.17900   
race2. Black | -0.85852 | 0.32331 | -2.655 | 0.00792 **
race3. Other | -0.70449 | 0.24624 | -2.861  | 0.00422 **

Signif. codes:  
0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

(Dispersion parameter for binomial family taken to be 1)

Null deviance: 603.05  on 486  degrees of freedom
Residual deviance: 570.53  on 482  degrees of freedom
AIC: 580.53

Number of Fisher Scoring iterations: 4

As we can see, after carrying out this multiple logistic regression, adjusting by age and race, we observe that both smoking and race are significantly (at 5% significance level) associated with low birthweight, whereas age is not. Thereby, we can infer that race is a confounder, while age is not. Despite the adjustment, the variance-covariance matrix of the effect estimates done by this model does not improve significantly, as shown by the vcov R command. However, since the values of this matrix are acceptable, this is not an issue :

```
vcov(model)[2,2]
```
> [1] 0.0399078

```
vcov(model1)[2,2]
```
>[1] 0.04635483


At this point, we have carried out a logistic regression in R and have identified potential confounders in the estimation of the effect of the variable smoker on the outcome (low birthweight). We thought of the possibility of doing this process multiple times, like in the “Z_X_Y_adjust_Z.R” script submitted by the professor on Moodle, but then we realised this would not be possible since, in that script, the professor defines the covariables of the linear model by the generation of random numbers, something that cannot be done in this situation since we are using a pre-existing dataset (clslowbwt) whose data cannot change. Therefore, we reckon we cannot do anything else with these logistic models, although we have proven that adjusting by different confounders in the model is totally possible in a logistic model.


# Bibliography

Matthew Blackwell. (2013). Observational Studies and Confounding. [Notes from Causal Inference couse, Harvard University] https://www.mattblackwell.org/files/teaching/s06-observational.pdf

Britannica, T. Information Architects of Encyclopaedia. (2022). thought. *Encyclopedia Britannica*. https://www.britannica.com/facts/thought

Foraita R., Spallek J., & Zeeb H. (2014). *Directed Acyclic Graphs. In: Ahrens W., Pigeot I. (eds) Handbook of Epidemiology*. Springer. https://doi.org/10.1007/978-0-387-09834-0_65

Pearl, Judea. (2009). Causal inference in statistics: An overview. *Statistics Surveys*, *3*, 96–146. doi:10.1214/09-SS057.

Pearl, J., Glymour, M., & Jewell., N. P. (2016). *Causal Inference in Statistics: A Primer, First Edition*. John Wiley & Sons Ltd.

Sjölander A. (2018). Estimation of causal effect measures with the R-package stdReg. *European journal of epidemiology*, *33*(9), 847–858. https://doi.org/10.1007/s10654-018-0375-y


