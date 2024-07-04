rm(list = ls(all = TRUE))
library(ggplot2)

# change to your owner folder
wdir = "/Users/Peter/Documents/"

# Parameter setting
# Initialize parameters and storage for history
mu1_init_0 <- 1
mu2_init_0 <- 5
pi_init_0 <- 0.3
sigma <- 1
iterations <- 500

mu1_true <- 2
mu2_true <- 10
pi1_true <- 0.6



# EM algorithm for a mixture of two normals with convergence plot
# Expectation step
e_step <- function(data, pi, mu1, mu2, sigma) {
  tau1 <- pi * dnorm(data, mean = mu1, sd = sigma)
  tau2 <- (1 - pi) * dnorm(data, mean = mu2, sd = sigma)
  gamma <- tau1 + tau2
  tau1<-tau1/gamma
  tau2<-tau2/gamma
  return(list(tau1=tau1, tau2=tau2, gamma=gamma))
}

# Maximization step
m_step <- function(data, tau1, tau2, gamma) {
  mu1 <- sum(data*tau1) / sum(tau1)
  mu2 <- sum(data*tau2) / sum(tau2)
  pi <- mean(tau1)
  return(list(mu1 = mu1, mu2 = mu2, pi = pi))
}

# Simulate data from a mixture of two normals
set.seed(0)
n = 250
x1 <- rnorm(pi1_true * n, mu1_true, sigma)
x2 <- rnorm((1 - pi1_true) * n, mu2_true, sigma)
data <- sample(c(x1, x2)) # Shuffling the data

history <- matrix(NA, nrow = iterations, ncol = 3)  # To store the parameter history

# Run EM algorithm and store the parameter history
mu1_init = mu1_init_0
mu2_init = mu2_init_0
pi_init = pi_init_0

for (i in 1:iterations) {
  para_e <- e_step(data, pi_init, mu1_init, mu2_init, sigma)
  params <- m_step(data, para_e$tau1, para_e$tau2, para_e$gamma)
  history[i, ] <- c(params$mu1, params$mu2, params$pi)
  mu1_init <- params$mu1
  mu2_init <- params$mu2
  pi_init <- params$pi
}

# Assuming the 'history' matrix is already populated with parameter values from the EM algorithm

# Open a PDF device to save the plots

png(paste0(wdir,paste0("Convergence_N_",iterations,"_Mu1_",mu1_init_0,"Mu2_",mu2_init_0,"Pro_",pi_init_0,".png")), 
    width = 1500, height = 600, bg = "transparent")
# Setup the layout for the plots
par(mfrow = c(1, 2), oma = c(0, 0, 0, 0), mar = c(4, 4, 2, 1), las = 1)

# Plot the convergence of means
plot(1:iterations, history[, 1], type = 'l', col = 'blue', ylim = range(c(history[, 1:2], mu1_true, mu2_true)), 
     xlab = 'Iteration', ylab = '', main = 'Convergence of Means',
     cex.axis = 1.6, cex.lab = 1.5, cex.main = 1.4, lwd = 2)
lines(1:iterations, history[, 2], col = 'red',lwd = 2)
abline(h = mu1_true, col = 'blue', lty = 2)
abline(h = mu2_true, col = 'red', lty = 2)


# Plot the convergence of proportion
plot(1:iterations, history[, 3], type = 'l', col = 'orange', ylim = range(c(history[, 3], pi1_true)),
     xlab = 'Iteration', ylab = '', main = 'Convergence of Proportion',
     cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.4, lwd = 2)
abline(h = pi1_true, col = 'orange', lty = 2)

# Close the device
dev.off()

