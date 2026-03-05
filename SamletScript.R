

# Packagess

library(tidyverse)
library(latex2exp)

# settings

options(scipen = 9999, digits = 5)
theme_set(theme_bw())
# options(pager = "cat")
# options(width = 10)
options(max.print = 10)

# data

mortality <- read.csv("mortality.csv", sep = ";", dec = ",")

mortality %>%
  ggplot(aes(x = Age, y = Intensity)) +
  geom_line()

mu_t <- approxfun(x = mortality$Age, y = mortality$Intensity, method = "linear", rule = 2)

# Parameters

r <- 0.02 
alpha <- 0.04
sigma <- 0.13
m <- 70
theta <- 0.12
K <- 4

# L(t) defined as a function of time t
L <- function(t) {
  600000 * exp(0.018 * (t - 30))
}

pi_t <- function(t) {
  ifelse(t < 55, 
         1,
         ifelse(t >= 55 & t <= m, 
                -(1 / (2 * (m - 55))) * t + (2 * m - 55) / (2 * (m - 55)), 
                0.5))
}

r_sim <- function(t) {
  pi_t(t)*alpha + (1-pi_t(t))*r
}

# Exercise (i)

thiele_ode <- function(t, V, r_t, mu_t) {
  # r_tilde and mu_tilde are functions of t from Table 2
  return((r_t(t) + mu_t(t)) * V - 1)
}

solve_v0_rk4 <- function(t_start, t_end, n, r_tilde_f, mu_tilde_f) {
  # We solve BACKWARD from t_end (111) to t_start (30)
  h <- (t_start - t_end) / n  # h will be negative
  t_vec <- seq(t_end, t_start, by = h)
  v_vec <- numeric(length(t_vec))
  
  v_vec[1] <- 0 # Boundary condition V(111) = 0
  
  for (i in 1:(length(t_vec) - 1)) {
    t <- t_vec[i]
    v <- v_vec[i]
    
    k1 <- h * thiele_ode(t, v, r_tilde_f, mu_tilde_f)
    k2 <- h * thiele_ode(t + h/2, v + k1/2, r_tilde_f, mu_tilde_f)
    k3 <- h * thiele_ode(t + h/2, v + k2/2, r_tilde_f, mu_tilde_f)
    k4 <- h * thiele_ode(t + h, v + k3, r_tilde_f, mu_tilde_f)
    
    v_vec[i+1] <- v + (k1 + 2*k2 + 2*k3 + k4) / 6
  }
  
  # Return as a function for easy lookup in the Euler scheme (part ii)
  return(approxfun(t_vec, v_vec))
}

V_func <- solve_v0_rk4(30, 111, 1000, r_sim, mu_t)

reserve_df <- data.frame(
  t = seq(30, 111, by = 0.01)
)
reserve_df$V <- V_func(reserve_df$t)

reserve_df

reserve_df %>%
  ggplot(aes(x = t, y = V)) +
  geom_line(size = 1.1) +
  scale_x_continuous(breaks = seq(30, 110, 10)) +
  scale_y_continuous(breaks = seq(0, 30, 5)) +
  xlab("t") +
  ylab(TeX("$\\tilde{V}_0(t)$"))

# Exercise (ii)

h <- 1
data.frame(t = seq(0, 111, 1)) %>%
  mutate(Z = rnorm(nrow(.)))

X_working <- function(t, h, z, x){
  (pi_t(t)*alpha + (1-pi_t(t))*r)*x*h
  + sigma*pi_t(t)*x*sqrt(h)*z
  + theta*L(t)*h
  - mu_t(t)*(K*L(t)-x)*h
}

simulate_account_path <- function(t_start = 30, t_end = 111, h = 1/10000, 
                                  X_init = 0, seed = NULL) {
  
  t_seq <- seq(t_start, t_end, by = h)
  N <- length(t_seq)
  
  X <- numeric(N)
  X[1] <- X_init 
  
  Z <- rnorm(N - 1)
  
  for (i in 1:(N - 1)) {
    t <- t_seq[i]
    x_curr <- X[i]
    
    diffusion <- sigma * pi_t(t) * x_curr
    
    if (t < m) {
      drift <- (pi_t(t) * alpha + (1 - pi_t(t)) * r) * x_curr + 
        theta * L(t) - 
        mu_t(t) * (K * L(t) - x_curr)
    } else {
      drift <- (pi_t(t) * alpha + (1 - pi_t(t)) * r) * x_curr - 
        (x_curr / V_func(t)) + 
        mu_t(t) * x_curr
    }
    
    X[i + 1] <- x_curr + drift * h + diffusion * sqrt(h) * Z[i]
  }
  
  return(tibble(Time = t_seq, Account_Value = X))
}

# simulated_path <- simulate_account_path(h = 1/10000)

# saveRDS(simulated_path, "SimPath.rds")
simulated_path <- readRDS("SimPath.rds")

simulated_path %>%
  mutate(Account_Value = Account_Value/10^6) %>%
  ggplot(aes(x = Time, y = Account_Value)) +
  geom_line() +
  scale_x_continuous(breaks = seq(0, 120, 10)) +
  geom_vline(xintercept = 70, linetype = "dashed") +
  xlab("t") +
  ylab("Account value (in millions)")

simulated_path$Account_Value

simulated_path %>%
  filter(Time == 70)
