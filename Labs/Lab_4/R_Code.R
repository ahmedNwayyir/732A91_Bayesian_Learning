#setwd("E:/1. Workshop/11. Bayesian Learning/1. Labs/Lab 4")
setwd("C:/Users/WizzCon/Desktop/Machine Learning/1. Workshop/11. Bayesian Learning/1. Labs/Lab 4")
################################################################################
## Time series models in Stan.
################################################################################

## @knitr a
## A

phi <- seq(-0.9,0.9, length.out = 9)

AR <- function(mu, sdSq, phi, T){
  X <- data.frame(matrix(NaN, nrow = T, ncol = length(phi)))
  X[1,] <- mu
  
  for(i in phi){
    for(t in 2:T){
      e <- rnorm(1, mean = 0, sd = sqrt(sdSq))
      X[t,which(phi == i)] <- mu + i * (X[t-1,which(phi == i)] - mu) + e
    }
  }
  return(data.frame(t = 1:T, X))
}

X <- AR(mu = 10, sdSq = 2, phi = phi, T = 200)
#colnames(X) <- c("t", "x1", "x2", "x3", "x4", "x5", "x6", "x7", "x8", "x9")

library(ggplot2)
library("latex2exp")
p1 <- ggplot(X) +
  geom_line(aes(x = X[,1], y = X[,2]), color = "#16264c", alpha = 0.7, size = 0.4) +
  labs(subtitle = TeX("$\\phi = -0.9$"), x = "t", y = "x") +
  theme_minimal()

p2 <- ggplot(X) +
  geom_line(aes(x = X[,1], y = X[,3]), color = "#16264c", alpha = 0.7, size = 0.4) +
  labs(subtitle = TeX("$\\phi = -0.675$"), x = "t", y = "x") +
  theme_minimal()

p3 <- ggplot(X) +
  geom_line(aes(x = X[,1], y = X[,4]), color = "#16264c", alpha = 0.7, size = 0.4) +
  labs(subtitle = TeX("$\\phi = -0.45$"), x = "t", y = "x") +
  theme_minimal()

p4 <- ggplot(X) +
  geom_line(aes(x = X[,1], y = X[,5]), color = "#16264c", alpha = 0.7, size = 0.4) +
  labs(subtitle = TeX("$\\phi = -0.225$"), x = "t", y = "x") +
  theme_minimal()

p5 <- ggplot(X) +
  geom_line(aes(x = X[,1], y = X[,6]), color = "#16264c", alpha = 0.7, size = 0.4) +
  labs(subtitle = TeX("$\\phi = 0$"), x = "t", y = "x") +
  theme_minimal()

p6 <- ggplot(X) +
  geom_line(aes(x = X[,1], y = X[,7]), color = "#16264c", alpha = 0.7, size = 0.4) +
  labs(subtitle = TeX("$\\phi = 0.225$"), x = "t", y = "x") +
  theme_minimal()

p7 <- ggplot(X) +
  geom_line(aes(x = X[,1], y = X[,8]), color = "#16264c", alpha = 0.7, size = 0.4) +
  labs(subtitle = TeX("$\\phi = 0.45$"), x = "t", y = "x") +
  theme_minimal()

p8 <- ggplot(X) +
  geom_line(aes(x = X[,1], y = X[,9]), color = "#16264c", alpha = 0.7, size = 0.4) +
  labs(subtitle = TeX("$\\phi = 0.675$"), x = "t", y = "x") +
  theme_minimal() 

p9 <- ggplot(X) +
  geom_line(aes(x = X[,1], y = X[,10]), color = "#16264c", alpha = 0.7, size = 0.4) +
  labs(subtitle = TeX("$\\phi = 0.9$"), x = "t", y = "x") +
  theme_minimal()


library(gridExtra)
grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, nrow = 3)


## @knitr b_i
## B_i

X <- AR(mu = 10, sdSq = 2, phi = c(0.3,0.95), T = 200)
suppressPackageStartupMessages(library(rstan))

# This step is optional, but it can result in compiled Stan programs that execute 
# much faster than they otherwise would. Simply paste the following into R once

# dotR <- file.path(Sys.getenv("HOME"), ".R")
# if (!file.exists(dotR)) dir.create(dotR)
# M <- file.path(dotR, ifelse(.Platform$OS.type == "windows", "Makevars.win", "Makevars"))
# if (!file.exists(M)) file.create(M)
# cat("\nCXX14FLAGS=-O3 -march=native -mtune=native",
#     if( grepl("^darwin", R.version$os)) "CXX14FLAGS += -arch x86_64 -ftemplate-depth-256" else 
#       if (.Platform$OS.type == "windows") "CXX11FLAGS=-O3 -march=corei7 -mtune=corei7" else
#         "CXX14FLAGS += -fPIC",
#     file = M, sep = "\n", append = TRUE)

# If you ever need to change anything with your C++ toolchain configuration, you can execute:
# M <- file.path(Sys.getenv("HOME"), ".R", ifelse(.Platform$OS.type == "windows", "Makevars.win", "Makevars"))
# file.edit(M)

# if you are using rstan locally on a multicore machine and have plenty of RAM 
# to estimate your model in parallel, at this point execute or we can set the number
# of cores in the stan() function
options(mc.cores = parallel::detectCores())

# allows you to automatically save a bare version of a compiled Stan program 
# to the hard disk so that it does not need to be recompiled (unless you change it).
rstan_options(auto_write = TRUE)


n_chains <- 4
n_cores  <- 4
warmups  <- 1000
iters    <- 2000
max_treedepth <- 10
accept_rate   <- 0.8

fit1 <- stan(file = "Other/Stan/AR1.stan", model_name = "AR1", data = list(N = dim(X)[1], x = X[,2]), 
             chains = n_chains, cores = n_cores, warmup = warmups, iter = iters, 
             control = list(max_treedepth = max_treedepth, adapt_delta = accept_rate))

fit2 <- stan(file = "Other/Stan/AR1.stan", model_name = "AR1", data = list(N = dim(X)[1], x = X[,3]), 
             chains = n_chains, cores = n_cores, warmup = warmups, iter = iters, 
             control = list(max_treedepth = max_treedepth, adapt_delta = accept_rate))


x_t <- extract(fit1, permuted = TRUE, inc_warmup = FALSE)
# The summaries for the parameters shown by the print method 
# are calculated using only post-warmup draws.

y_t <- extract(fit2, permuted = TRUE, inc_warmup = FALSE)

plot(fit1, ci_level = 0.95) +
  labs(subtitle = TeX("$\\x_t$ 95% Credible Interval"), x = "", y = "") +
  theme_bw()


plot(fit2, ci_level = 0.95) +
  labs(subtitle = TeX("$\\y_t$ 95% Credible Interval"), x = "", y = "") +
  theme_bw()

library(knitr)
kable(summary(fit1, pars = c("mu", "phi", "sdSq"), probs = c(0.025, 0.975))$summary[,c(1,3,4,5)], 
             caption = "x 95% Credible Interval")
kable(summary(fit2, pars = c("mu", "phi", "sdSq"), probs = c(0.025, 0.975))$summary[,c(1,3,4,5)], 
             caption = "y 95% Credible Interval")

## @knitr b_ii_1
## B_ii

getPlotData <- function(model, n_chains, iters, warmups){
  post_warmup <- iters - warmups
  
  chains <- NULL
  for(chain in 1:n_chains){
    Cumsum <- apply(As.mcmc.list(model)[[chain]], 2, function(x)cumsum(x)/ seq(1,post_warmup))
    colnames(Cumsum) <- c("mu_cumsum", "phi_cumsum", "sdSq_cumsum", "lp_cumsum")
    single_chain <- cbind(xGrid = 1:post_warmup, As.mcmc.list(model)[[chain]], Cumsum,
                             div = get_sampler_params(model, inc_warmup = FALSE)[[chain]][,5])
    chains <- rbind(chains, single_chain)
  }

  return(as.data.frame(chains))
}

x_plotData <- getPlotData(fit1, n_chains, iters, warmups)
y_plotData <- getPlotData(fit2, n_chains, iters, warmups)
x_divDraws   <- data.frame(mu = x_plotData$mu[which(x_plotData$div == 1)], 
                           phi = x_plotData$phi[which(x_plotData$div == 1)])
y_divDraws   <- data.frame(mu = y_plotData$mu[which(y_plotData$div == 1)], 
                           phi = y_plotData$phi[which(y_plotData$div == 1)])


##### Joint Density Plots #####
ggplot(x_plotData)+
  geom_point(aes(x = mu, y = phi), color = "#7475FD", alpha = 0.5, size = 0.6) +
  geom_point(data = x_divDraws, aes(x = mu, y = phi), color = "#F64769", alpha = 0.3, size = 0.7) +
  labs(subtitle = TeX("Posterior Distribution of ($\\mu_1$, $\\phi_1$)"),
       x = TeX("$\\mu_1$"), y = TeX("$\\phi_1$")) +
  theme_bw()

ggplot(y_plotData)+
  geom_point(aes(x = mu, y = phi), colour = "#7475FD", alpha = 0.5, size = 0.6) +
  geom_point(data = y_divDraws, aes(x = mu, y = phi), colour = "red", alpha = 0.3, size = 0.7) +
  labs(subtitle = TeX("Posterior Distribution of ($\\mu_2$, $\\phi_2$)"),
       x = TeX("$\\mu_2$"), y = TeX("$\\phi_2$")) +
  theme_bw()


## @knitr b_ii_2
##### x Convergance Plots #####
gp1 <- ggplot(x_plotData) +
  geom_line(aes(x = xGrid, y = mu), color = "red", alpha = 0.2, size = 0.2) +
  geom_point(aes(x = xGrid, y = mu), color = "red", alpha = 0.2, size = 0.5) +
  geom_line(aes(x = xGrid, y = mu_cumsum), color = "black", alpha = 0.7, size = 0.2) +
  labs(subtitle = TeX("Convergence of $\\mu_x$"), x = "", y = TeX('$\\mu_x$')) +
  theme_minimal() +
  theme(legend.position="bottom", legend.title = element_blank())

gp2 <- ggplot(x_plotData) +
  geom_line(aes(x = xGrid, y = phi), color = "#FCFA2922", alpha = 0.2, size = 0.2) +
  geom_point(aes(x = xGrid, y = phi), color = "#FCFA2922", alpha = 0.2, size = 0.5) +
  geom_line(aes(x = xGrid, y = phi_cumsum), color = "black", alpha = 0.7, size = 0.2) +
  labs(subtitle = TeX("Convergence of $\\phi_x$"), x = "", y = TeX('$\\phi_x$')) +
  theme_minimal() +
  theme(legend.position="bottom", legend.title = element_blank())

gp3 <- ggplot(x_plotData) +
  geom_line(aes(x = xGrid, y = sdSq), color = "#7475FD", alpha = 0.2, size = 0.2) +
  geom_point(aes(x = xGrid, y = sdSq), color = "#7475FD", alpha = 0.2, size = 0.5) +
  geom_line(aes(x = xGrid, y = sdSq_cumsum), color = "black", alpha = 0.7, size = 0.2) +
  labs(subtitle = TeX("Convergence of $\\sigma_x^2$"), x = "", y = TeX('$\\sigma_x^2$')) +
  theme_minimal() +
  theme(legend.position="bottom", legend.title = element_blank())

gp4 <- ggplot(y_plotData) +
  geom_line(aes(x = xGrid, y = mu), color = "red", alpha = 0.2, size = 0.2) +
  geom_point(aes(x = xGrid, y = mu), color = "red", alpha = 0.2, size = 0.5) +
  geom_line(aes(x = xGrid, y = mu_cumsum), color = "black", alpha = 0.7, size = 0.2) +
  labs(subtitle = TeX("Convergence of $\\mu_y$"), x = "", y = TeX('$\\mu_y$')) +
  theme_minimal() +
  theme(legend.position="bottom", legend.title = element_blank())

gp5 <- ggplot(y_plotData) +
  geom_line(aes(x = xGrid, y = phi), color = "#FCFA2922", alpha = 0.2, size = 0.2) +
  geom_point(aes(x = xGrid, y = phi), color = "#FCFA2922", alpha = 0.2, size = 0.5) +
  geom_line(aes(x = xGrid, y = phi_cumsum), color = "black", alpha = 0.7, size = 0.2) +
  labs(subtitle = TeX("Convergence of $\\phi_y$"), x = "", y = TeX('$\\phi_y$')) +
  theme_minimal() +
  theme(legend.position="bottom", legend.title = element_blank())

gp6 <- ggplot(y_plotData) +
  geom_line(aes(x = xGrid, y = sdSq), color = "#7475FD", alpha = 0.2, size = 0.2) +
  geom_point(aes(x = xGrid, y = sdSq), color = "#7475FD", alpha = 0.2, size = 0.5) +
  geom_line(aes(x = xGrid, y = sdSq_cumsum), color = "black", alpha = 0.7, size = 0.2) +
  labs(subtitle = TeX("Convergence of $\\sigma_y^2$"), x = "", y = TeX('$\\sigma_y^2$')) +
  theme_minimal() +
  theme(legend.position="bottom", legend.title = element_blank())

grid.arrange(gp1, gp4, gp2, gp5, gp3, gp6, ncol = 2)


## @knitr b_ii_3
# The following code is taken from this blog with slight adjustment
# https://www.weirdfishes.blog/blog/fitting-bayesian-models-with-stan-and-r/
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(dplyr))
diagnostics_1 <- get_sampler_params(fit1) %>% 
  set_names(1:n_chains) %>% 
  map_df(as_tibble,.id = 'chain') %>% 
  group_by(chain) %>% 
  mutate(iteration = 1:length(chain)) %>% 
  mutate(warmup = iteration <= warmups) 

diagnostics_1 %>% 
  group_by(warmup, chain) %>% 
  summarise(percent_divergent = mean(divergent__ >0)) %>% 
  ggplot() +
  geom_col(aes(chain, percent_divergent, fill = warmup), 
           position = 'dodge', color = 'black') + 
  scale_y_continuous(labels = scales::percent, name = "% Divergent Runs") +
  scale_fill_manual(values = c("#16264c", "#ACBEE7")) +
  labs(subtitle = "% of Divergent Transitions for x_t") +
  theme_bw() +
  theme(legend.position="bottom")

diagnostics_2 <- get_sampler_params(fit2) %>% 
  set_names(1:n_chains) %>% 
  map_df(as_tibble,.id = 'chain') %>% 
  group_by(chain) %>% 
  mutate(iteration = 1:length(chain)) %>% 
  mutate(warmup = iteration <= warmups) 

diagnostics_2 %>% 
  group_by(warmup, chain) %>% 
  summarise(percent_divergent = mean(divergent__ >0)) %>% 
  ggplot() +
  geom_col(aes(chain, percent_divergent, fill = warmup), 
           position = 'dodge', color = 'black') + 
  scale_y_continuous(labels = scales::percent, name = "% Divergent Runs") +
  scale_fill_manual(values = c("#16264c", "#ACBEE7")) +
  labs(subtitle = "% of Divergent Transitions for y_t") +
  theme_bw() +
  theme(legend.position="bottom")



## @knitr c
## C

campy <- read.table("Other/Data/campy.dat", header=TRUE)

n_chains <- 4
n_cores  <- 4
warmups  <- 1000
iters    <- 4000
max_treedepth <- 10
accept_rate   <- 0.8

ct_fit <- stan(file = "Other/Stan/AR2.stan", model_name = "AR2", data = list(N = dim(campy)[1], y = campy), 
               chains = n_chains, cores = n_cores, warmup = warmups, iter = iters, 
               control = list(max_treedepth = max_treedepth, adapt_delta = accept_rate))



ct <- extract(ct_fit)
ct_summ <- summary(ct_fit)$summary

mean_ci <- exp(summary(ct_fit, probs = c(0.025, 0.975))$summary[c(4:143),c(1,4,5)])
ctPlotData <- data.frame(Time = 1:nrow(campy), 
                         Data  = campy$c, 
                         Mean  = mean_ci[,1],
                         Lower = mean_ci[,2],
                         Upper = mean_ci[,3])

ggplot(ctPlotData, aes(x = Time)) +
  geom_point(aes(y = Data, col = "Data"), alpha = 0.6) + 
  geom_line(aes(y = Mean, col = "Mean"), size = 0.5) +
  geom_line(aes(y = Lower, col = "95% Bounds"), lty = 2, size = 0.4, alpha = 0.8) +
  geom_line(aes(y = Upper), color = "blue", lty = 2, size = 0.4, alpha = 0.8) +
  scale_color_manual(breaks = c("Data", "Mean", "95% Bounds"), 
                     values = c("black", "red", "blue")) +
  labs(y = TeX("$\\c_t$")) +
  theme_bw() +
  theme(legend.position="bottom", legend.title = element_blank())


## @knitr d
## D

ct2_fit <- stan(file = "Other/Stan/AR3.stan", model_name = "AR3", data = list(N = dim(campy)[1], y = campy), 
               chains = n_chains, cores = n_cores, warmup = warmups, iter = iters, 
               control = list(max_treedepth = max_treedepth, adapt_delta = accept_rate))



ct2 <- extract(ct2_fit)
ct2_summ <- summary(ct2_fit)$summary

mean_ci_2 <- exp(summary(ct2_fit, probs = c(0.025, 0.975))$summary[c(4:143),c(1,4,5)])
ct2PlotData <- data.frame(Time = 1:nrow(campy), 
                          Data  = campy$c, 
                          Mean  = mean_ci_2[,1],
                          Lower = mean_ci_2[,2],
                          Upper = mean_ci_2[,3])

ggplot(ct2PlotData, aes(x = Time)) +
  geom_point(aes(y = Data, col = "Data"), alpha = 0.6) + 
  geom_line(aes(y = Mean, col = "Mean"), size = 0.5) +
  geom_line(aes(y = Lower, col = "95% Bounds"), lty = 2, size = 0.4, alpha = 0.8) +
  geom_line(aes(y = Upper), color = "blue", lty = 2, size = 0.4, alpha = 0.8) +
  scale_color_manual(breaks = c("Data", "Mean", "95% Bounds"), 
                     values = c("black", "red", "blue")) +
  labs(y = TeX("$\\c_t$")) +
  theme_bw() +
  theme(legend.position="bottom", legend.title = element_blank())