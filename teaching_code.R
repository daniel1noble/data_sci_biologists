

# We want really large sample sizes so that the analytical solutions (for SE, for example) converge on the sampled data
n = 10000

# Simulate two groups, one with mean = 5 and otehr mean = 10
x <- rnorm(n, 5, 1)
y <- rnorm(n, 10, 1)

# Calculate their standard errors
se_x <- sd(x) / sqrt(length(x))
se_y <- sd(y) / sqrt(length(y))

# Or sampling variance
sd(y)^2 / length(y) == se_y^2

# Now, calculate the difference
diff <- y - x

# Calculate the SE of the difference. Remember, this is the analytical solution, it will always be exact
diff_se <- sd(diff) / sqrt(length(diff))

# Calculate SE of diffs from individual group SE's. This should pretty much match diff_se
diff_se_2 <- sqrt(se_x^2 + se_y^2)

# Just playing with ocode for lnRR 
## Log response ratio 
   lnrr <- log(y / x)              # Log RR for y and x
   lnRR <- log(mean(y) / mean(x))  # Analytical solution using the mean of y and x

# Log response ratio sampling variance
 sv_lrr <- (sd(lnrr) / sqrt(length(lnrr)))^2 # Sampling variance using SE
  vlnRR <- sd(y)^2 / (mean(y)^2*length(y)) + sd(x)^2 / (mean(x)^2*length(x)) # Analytical solution for sampling variance