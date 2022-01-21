# Programming-assignment: Illustrating problems in the adjustment of covariates when doing causal inference
Authors: Patricia Fernandez Moreno, Sofía González Matatoros, Víctor Manuel López Molina and Álvaro Pita Ramírez

# 1. Introduction to causal inference
According to the Encyclopedia Britannica (2022), induction is a means of reasoning from a particulars to generals and one example of this is causal inference. The objective of this process is to infer the actual effect of a particular element within a more complex system, considering the global outcome under changing conditions (Pearl & Judea, 2009). At this point, it is important to clarify that association does not imply causation, in other words, two elements can be related because they belong to the same interacting network, but its relationship cannot be due to a cause-effect process.

Back to the initial topic, when we want to make an interpretation of the statistical results in this type of scenarios, a relevant aspect is to choose on which covariates we are going to adjust in the analysis. Therefore, the aim of this repport is to explain which are the consequences of doing an incorrect analysis, because we are concerned that the process of causal inference is highly used in medical research.

# 2. Explaining the appropriate analysis

To illustrate how to perform the correct analysis, we are going to examinate the following scenario.

#incluir la imagen

```

#incluir código ajustado (Víctor)
dag1 <- dagitty("dag {
X -> Z
X -> Y
Y -> W
e_z -> Z
e_y -> Y
e_w -> W
}")

coordinates(dag1) <- list(x = c(Z = 1, X = 2, Y = 3, e_x = 1.75, e_y = 2.75),
                          y = c(Z = 0, X = 0, Y = 0, e_x = -.2, e_y = -.2))

## plot(dag1)
drawdag(dag1) #se dibuja el DAG
drawdag(dag1, xlim = c(0.5, 3.5), ylim = c(-2, 1)) #se dibuja el DAG pero más pequeño

```
Considering this DAG, we can make these observations
- Z and Y are associated throw X, which is a cofounder, so conditioning on X would render them independent.
- X and W are connected throw Y, but conditioning on Y will render them independent.
- Z and W are associated throw X and Y.

In conclusion, to analyze the data related to this plot we should not adjust by any covariate. 

In order to have a reference, we are going to calculate the summary and the variance of the regression models applied to this DAG.

```
summary(reg1) # Z ~ X
vcov(reg1)[2,2] 

summary(reg2) # Y ~ X
vcov(reg2)[2,2] 

summary(reg3) # W ~ Y
vcov(reg3)[2,2]

summary(reg4) # W ~ X + Z
vcov(reg4)[2,2]

```


# 3. Incorrect analysis

# 4. Comparisons and conclusions

# Bibliography
Britannica, T. Information Architects of Encyclopaedia (2022). thought. Encyclopedia Britannica. https://www.britannica.com/facts/thought
Pearl, Judea (2009). "Causal inference in statistics: An overview". Statistics Surveys. 3: 96–146. doi:10.1214/09-SS057.
