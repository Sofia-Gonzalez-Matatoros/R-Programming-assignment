
##LINEAR REGRESSION

#DAG of our model

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

##Finding confounders in a multiple linear regression.
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

dist_multiple()

##Finding confounders in a simple linear regression.
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

dist_simple()

##After creating these models, we changed some parameters of the
#function with the objective of trying to improve the adjustment by a certain
#confounder when this adjustment was not favourable in the initial situation.

#Specifically, we changed X's standard deviation and its coefficient in the
#multiple linear regression model, showing a better result when adjusting by W.
dist_multiple(sd_x = 5, b_zx = -4)

#The same result can be seen as for the simple linear regression.
dist_simple(sd_x = 5, b_zx = -4)

##LOGISTIC REGRESSION
#We performed a logistic regression by the use of the clslowbwt dataset.

#First, we had to install the AF package and load this dataset.
#install.packages("AF")
library(AF)
data(clslowbwt)

#Simple logistic regression model (without adjusting by any confounder)
model<-glm(formula=lbw~smoker,family="binomial",data=clslowbwt)
summary(model)

#Multiple logistic regression model (adjusting by race and age)
model1 <- glm(formula=lbw~smoker+age+race,family="binomial",data=clslowbwt)
summary(model1)

#Analysis of the variance in both models
vcov(model)[2,2]
vcov(model1)[2,2]



