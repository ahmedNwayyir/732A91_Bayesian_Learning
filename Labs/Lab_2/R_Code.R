setwd("C:/Users/WizzCon/Desktop/Machine_Learning/1_Workshop/9_Bayesian_Learning/1_Labs/Lab_2")

## @knitr 1_a1
################################################################################
## 1. Linear and polynomial regression
################################################################################

## A

link_temp <- read.table("Other/Data/TempLinkoping.txt", header=TRUE)
attach(link_temp)

mu_0  <- c(-10,100,-100)
om_0  <- diag(rep(0.01,3))
nu_0  <- 4
var_0 <- 1

inv_chi <- function(df, Var){
  chi <- rchisq(1, df)
  Var <- (df * Var) / chi
  return(Var)
}
model_var_0 <- inv_chi(df=nu_0, Var=var_0)


library(mvtnorm)
nPriors <- 100
X    <- cbind(rep(1,length(time)), time, time^2)
y    <- temp

prior_generator <- function(n, mu, om, model_var){
  yHat <- matrix(0, nrow=length(time), ncol=n)
  Xy   <- NULL
  for(i in 1:n){
    Beta      <- rmvnorm(n=1, mean=mu, sigma=model_var * solve(om))
    yHat[,i]  <- X %*% matrix(Beta)
    temporary <- data.frame(x=1:length(time), y=yHat[,i], col=rep(i:i, length(time)))
    Xy        <- rbind(Xy, temporary)
  }
  return(data.frame(Xy))
}
prior_data <- prior_generator(n=nPriors, mu=mu_0, om=om_0, model_var=model_var_0)

library(ggplot2)
ggplot(prior_data) +  
  geom_line(aes(x=x, y=y, group=col, color = "Priors"), alpha=0.4) +
  labs(subtitle="Priors from from the first set of hyperparameters", x = "Time", y = "Temperature") +
  scale_color_manual(values = c("red")) +
  theme_bw() +
  theme(legend.position="bottom", legend.title = element_blank()) 

## @knitr 1_a2
# The intercept mean (starting day which is in January) changed to reflect
# our belief of what might be the mean temperature in the first 2 months
# Changed B_1 & B_2 to 98 & -96 respectively to make the middle hump higher "temperatures in summer between 10-30
# and the end of the year similar to the begining of the year
mu_0  <- c(-3,98,-96)    
# The om_0 increased to 0.15 to improve the accuracy otherwise we could reduce the variance hyper-parameter
om_0  <- diag(c(0.15,0.15,0.15))
# Increased nu to 100 to give us more control of the outcome of model variance "more certainty"
nu_0  <- 100
var_0 <- 1
model_var_0 <- inv_chi(df=nu_0, Var=var_0)

prior_data <- prior_generator(n=nPriors, mu=mu_0, om=om_0, model_var=model_var_0)

ggplot(prior_data) +  
  geom_line(aes(x=x, y=y, group=col, color = "Priors"), alpha=0.4) +
  labs(subtitle="Priors from from the second set of hyperparameters", x = "Time", y = "Temperature") +
  scale_color_manual(values = c("red")) +
  theme_bw() +
  theme(legend.position="bottom", legend.title = element_blank()) 

## @knitr 1_a3
ggplot(prior_data) +  
  geom_point(data=link_temp, aes(x=1:length(time), y=temp, color="Temperature Data"), alpha=0.4) +
  labs(subtitle="Linkoping Temperature Data", x = "Time", y = "Temperature") +
  scale_color_manual(values = c("red")) +
  theme_bw() +
  theme(legend.position="bottom", legend.title = element_blank()) 


## @knitr 1_b1
## B

n    <- length(time)
B_h  <- as.numeric(lm(temp ~ time + I(time^2), data=link_temp)$coefficients)
mu_n <- as.numeric(solve(t(X) %*% X + om_0) %*% (t(X) %*% X %*% B_h + om_0 %*% mu_0))
om_n <- matrix(t(X) %*% X + om_0, ncol=3)
nu_n <- nu_0 + n
nu_n_var_0 <- nu_0 * var_0 + (t(y) %*% y + (t(mu_0) %*% om_0 %*% mu_0) - (t(mu_n) %*% om_n %*% mu_n))
var_n <- as.numeric(nu_n_var_0/nu_n)

nPosters <- 1000
posterPmeters <- NULL
for (i in 1:nPosters) {
  model_var_n   <- inv_chi(df=nu_n, Var=var_n)
  B_n           <- rmvnorm(n=1, mean=mu_n, sigma=model_var_n * solve(om_n))
  temporary     <- data.frame(Var=model_var_n, B_0=B_n[[1]], B_1=B_n[[2]], B_2=B_n[[3]])
  posterPmeters <- rbind(posterPmeters, temporary)
}

p1 <- ggplot(posterPmeters) +
  geom_histogram(aes(x=Var, y=..density..), bins=20, fill="#ffffffff", colour="black", size=0.1) +
  labs(subtitle="Model Posterior Variance", y="Density", x="Variance") +
  theme_minimal()
p2 <- ggplot(posterPmeters) +
  geom_histogram(aes(x=B_0, y=..density..), bins=20, fill="#ffffffff", colour="black", size=0.1) +
  labs(subtitle="Beta_0", y="Density", x="Beta_0") +
  theme_minimal()
p3 <- ggplot(posterPmeters) +
  geom_histogram(aes(x=B_1, y=..density..), bins=20, fill="#ffffffff", colour="black", size=0.1) +
  labs(subtitle="Beta_1", y="Density", x="Beta_1") +
  theme_minimal()
p4 <- ggplot(posterPmeters) +
  geom_histogram(aes(x=B_2, y=..density..), bins=20, fill="#ffffffff", colour="black", size=0.1) +
  labs(subtitle="Beta_2", y="Density", x="Beta_2") +
  theme_minimal()

library(gridExtra)
grid.arrange(p1, p2, p3, p4, nrow = 2)

## @knitr 1_b2
Beta <- matrix(0, nrow = nPosters, ncol = 3)
poster_generator <- function(n, mu, om, model_var){
  yHat <- matrix(0, nrow=length(time), ncol=n)
  for(i in 1:n){
    Beta[i,] <- rmvnorm(n=1, mean=mu, sigma=model_var * solve(om))
    yHat[,i] <- X %*% matrix(Beta[i,])
  }
  return(list(Beta, yHat))
}

poster_model  <- poster_generator(n=nPosters, mu=mu_n, om=om_n, model_var=model_var_n)[[2]]

poster_mean   <- apply(poster_model, 1, mean)
poster_median <- apply(poster_model, 1, median)
lower_bound   <- apply(poster_model, 1, quantile, probs= c(0.025, 0.975))[1,]
upper_bound   <- apply(poster_model, 1, quantile, probs= c(0.025, 0.975))[2,]
 
poster_data <- data.frame(time=1:length(time), temp, poster_mean, poster_median, lower_bound, upper_bound)
ggplot(poster_data) +
  geom_point(aes(x=time, y=temp), color="red", alpha=0.4) +
  geom_ribbon(aes(x=time, ymin=lower_bound, ymax=upper_bound, color="Credible Interval"), alpha=0.05) +
  geom_line(aes(x=time, y=poster_median, color="Posterior Median"), size = 1, linetype = 2) +
  labs(subtitle="Posterior Median & Credible Interval", y="Temperature", x="Time") +
  scale_color_manual(breaks=c("Credible Interval", "Posterior Median"), 
                     values=c("#AF8969", "#856347")) +
  theme_bw() +
  theme(legend.position="bottom", legend.title=element_blank()) 


## @knitr 1_c
## C
Betas   <- poster_generator(n=nPosters, mu=mu_n, om=om_n, model_var=model_var_n)[[1]]
x_tilde <- -Betas[,2] / (2 * Betas[,3])
dates   <- as.Date(365*x_tilde, origin="2018-01-01") 

ggplot(as.data.frame(dates)) +
  geom_histogram(aes(x=dates, y=..density..), bins=20, color="black", fill="#ffffffff", size=0.2) +
  geom_density(aes(x=dates, y=..density..), color = "red", size = 0.7) +
  labs(subtitle="Posterior Distribution of x_tilde", y="Density", x="Date") +
  theme_bw()


## @knitr 1_d
## D
nPriors <- 100
X    <- cbind(rep(1,length(time)), time, time^2, time^3, time^4, time^5, time^6, time^7)
y    <- temp

mu_0  <- c(5, 98, -98, 0, 0, 0, 0, 0)    
om_0  <- diag(0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01)
nu_0  <- 100
var_0 <- 1
model_var_0 <- inv_chi(df=nu_0, Var=var_0)

prior_generator <- function(n, mu, om, model_var){
  yHat <- matrix(0, nrow=length(time), ncol=n)
  Xy   <- NULL
  for(i in 1:n){
    Beta      <- rmvnorm(n=1, mean=mu, sigma=model_var * solve(om))
    yHat[,i]  <- X %*% matrix(Beta)
    temporary <- data.frame(x=1:length(time), y=yHat[,i], col=rep(i:i, length(time)))
    Xy        <- rbind(Xy, temporary)
  }
  return(data.frame(Xy))
}
prior_data <- prior_generator(n=nPriors, mu=mu_0, om=om_0, model_var=model_var_0)
prior_data <- cbind(time=1:length(time), temp, prior_data)

ggplot(prior_data) +  
  geom_line(aes(x=x, y=y, group=col, color = "Priors"), alpha=0.4) +
  geom_point(aes(x=time, y=temp), color="red", alpha=0.4) +
  labs(subtitle="Priors from from the second set of hyperparameters", x = "Time", y = "Temperature") +
  scale_color_manual(values = c("red")) +
  theme_bw() +
  theme(legend.position="bottom", legend.title = element_blank()) 



n    <- length(time)
B_h  <- as.numeric(lm(temp ~ time + I(time^2) + I(time^3) + I(time^4) + I(time^5) + I(time^6) + I(time^7), data=link_temp)$coefficients)
mu_n <- as.numeric(solve(t(X) %*% X + om_0) %*% (t(X) %*% X %*% B_h + om_0 %*% mu_0))
om_n <- matrix(t(X) %*% X + om_0, ncol=8)
nu_n <- nu_0 + n
nu_n_var_0 <- nu_0 * var_0 + (t(Y) %*% Y + (t(mu_0) %*% om_0 %*% mu_0) - (t(mu_n) %*% om_n %*% mu_n))
var_n <- as.numeric(nu_n_var_0/nu_n)

Beta <- matrix(0, nrow = nPosters, ncol = 8)
poster_generator <- function(n, mu, om, model_var){
  yHat <- matrix(0, nrow=length(time), ncol=n)
  for(i in 1:n){
    Beta[i,] <- rmvnorm(n=1, mean=mu, sigma=model_var * solve(om))
    yHat[,i] <- X %*% matrix(Beta[i,])
  }
  return(list(Beta, yHat))
}

poster_model  <- poster_generator(n=nPosters, mu=mu_n, om=om_n, model_var=model_var_n)[[2]]

poster_mean   <- apply(poster_model, 1, mean)
poster_median <- apply(poster_model, 1, median)
lower_bound   <- apply(poster_model, 1, quantile, probs= c(0.025, 0.975))[1,]
upper_bound   <- apply(poster_model, 1, quantile, probs= c(0.025, 0.975))[2,]

poster_data <- data.frame(time=1:length(time), temp, poster_mean, poster_median, lower_bound, upper_bound)
ggplot(poster_data) +
  geom_point(aes(x=time, y=temp), color="red", alpha=0.4) +
  geom_ribbon(aes(x=time, ymin=lower_bound, ymax=upper_bound, color="Credible Interval"), alpha=0.05) +
  geom_line(aes(x=time, y=poster_median, color="Posterior Median"), size = 1, linetype = 2) +
  labs(subtitle="Posterior Median & Credible Interval", y="Temperature", x="Time") +
  scale_color_manual(breaks=c("Credible Interval", "Posterior Median"), 
                     values=c("#AF8969", "#856347")) +
  theme_bw() +
  theme(legend.position="bottom", legend.title=element_blank()) 


## @knitr 2_a
################################################################################
## 1. Linear and polynomial regression
################################################################################

women <- read.table("Other/Data/WomenWork.dat", header=TRUE)
tau   <- 10   

y <- as.vector(women[,1])
X <- as.matrix(women[,2:9])
nPara <- dim(X)[2]

mu    <- as.vector(rep(0, nPara))
Sigma <- tau^2*diag(nPara)

LogPostLogistic <- function(betaVect, y, X, mu, Sigma){
  nPara   <- length(betaVect)
  linPred <- X %*% betaVect
  
  logLik <- sum(linPred*y -log(1 + exp(linPred)))
  if (abs(logLik) == Inf) logLik = -20000 # Likelihood is not finite, stear the optimizer away from here!
  
  logPrior <- dmvnorm(betaVect, matrix(0, nPara, 1), Sigma, log=TRUE);
  return(logLik + logPrior)
}

initVal <- as.vector(rep(0,dim(X)[2]))
OptimResults <-optim(initVal,LogPostLogistic,gr=NULL,y,X,mu,Sigma,method=c("BFGS"),
                     control=list(fnscale=-1),hessian=TRUE)

postMode <- OptimResults$par
postCov  <- -solve(OptimResults$hessian)
names(postMode) <- colnames(X) 
colnames(postCov) <- rownames(postCov) <- colnames(X) 

glmModel <- glm(Work ~ 0 + ., data = women, family = binomial)

knitr::kable(cbind("Maximum likelihood" = glmModel$coefficients, "Posterior Mode" = postMode))
knitr::kable(postCov, digits = 4)

approxPostStd <- sqrt(diag(postCov))


betaValues <- seq(postMode[7] - 3*approxPostStd[7], postMode[7] + 3*approxPostStd[7], length = 10)
betaSample <- rnorm(1000, postMode[7], approxPostStd[7])

q1 <- quantile(betaSample,.025)
q2 <- quantile(betaSample,.975)
dens <- density(betaSample)
dens_df <- data.frame(x = dens$x, y = dens$y)

ggplot(as.data.frame(betaSample)) +
  geom_histogram(aes(x=betaSample, y=..density..), bins=30, color="black", fill="#ffffffff", size=0.1) +
  geom_density(aes(x=betaSample, y=..density..), color = "red", size = 0.7) +
  geom_area(data = subset(dens_df, x >= q1 & x <= q2), 
            aes(x=x,y=y), fill = 'red', alpha = 0.3) +
  labs(subtitle="Posterior Distribution of NSmallChild", y="Density", x="Beta") +
  theme_bw() +
  theme(legend.position="bottom", legend.title = element_blank()) 

cat("Equal Tail Interval:", q1, "-", q2)


## @knitr 2_b
X_pred <- matrix(c(1, 10, 8, 10, (10/10)^2, 40, 1, 1))

simulator <- function(n, X, mu, Sigma) {
  beta     <- rmvnorm(n, mean = mu, sigma = Sigma)
  logistic <- exp(beta %*% X) / (1 + exp(beta %*% X))
  #return(plogis(beta %*% X))
  return(logistic)
}

pred_dist <- simulator(n = 1000, X = X_pred, mu = postMode, Sigma = postCov)
pred_dens <- density(pred_dist)
pred_dens_df <- data.frame(pred_dens$x, pred_dens$y)

ggplot(as.data.frame(pred_dist)) +
  geom_histogram(aes(x = pred_dist, y=..density..),bins = 30, color = "black", fill = "#ffffffff", size = 0.1) +
  geom_line(data = pred_dens_df, aes(x = pred_dens.x, y = pred_dens.y), color = "red") +
  labs(subtitle="Posterior Predictive Distribution", y="Density", x="Response Variable") +
  theme_bw() +
  theme(legend.position="bottom", legend.title = element_blank()) 


## @knitr 2_c
X_pred <- matrix(c(1, 10, 8, 10, (10/10)^2, 40, 1, 1))

binom_sim  <- function(n, X, mu, Sigma) {
  logistic <- vector(length = dim(X)[1])
  binom    <- vector(length = n)
  for(i in 1:n){
    beta     <- rmvnorm(i, mean = mu, sigma = Sigma)
    logistic <- exp(beta %*% X) / (1 + exp(beta %*% X))
    trials   <- rbinom(n = 10, size = 1, prob = logistic)
    binom[i] <- sum(trials) / dim(X)[1]
  }
  return(binom)
}

binom_dist <- binom_sim(n = 1000, X = X_pred, mu = postMode, Sigma = postCov)
binom_dens <- density(binom_dist)
binom_dens_df <- data.frame(binom_dens$x, binom_dens$y)

ggplot(as.data.frame(binom_dist)) +
  geom_histogram(aes(x = binom_dist), bins = 30, color = "#113B69", fill = "#3486DF", size = 0.1) +
  labs(subtitle="Predictive distribution from 1000 trials of 10 women", y="Number of successful trials", x="Probability of working women") +
  theme_bw() +
  theme(legend.position="bottom", legend.title = element_blank()) 
