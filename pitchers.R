# This file contains the R code to conduct the experiments described in our report and re-construct the plots
# The code is divided section wise, with the sections same in the report.
################## Section 1: Introduction ####################

#################################################


################## Section 2: Introduction to OLLLTN ####################

set.seed(511)
library(pracma) #For logit function
G <- function(mu, sigma, y){ #Function that returns eita = G(y) for given mu and sigma
  q <- (logit(y)-logit(mu))/sigma 
  eita <- pnorm(q, 0.0, 1.0)
  return (eita)
}

Pdf <- function(mu, sigma, v, y){ #Returns the value of pdf of OLLLTN at y for given mu, sigma and nu
  p1 <- v/(y*(1-y)*(2*pi*(sigma)^2)^0.5)
  p2 <- exp(-((logit(y)-logit(mu))^2/(2*sigma^2)))
  eita <- G(mu, sigma, y)
  p3 <- (eita*(1-eita))^(v-1)
  p4 <- (eita^v + (1-eita)^v)^(-2)
  return (p1*p2*p3*p4)
}

Quantile <-function(u, mu, sigma, v){ #Returns the value of quantile of OLLLTN at u for given mu,sigma,nu
  u1 <- u^(1/v)
  u2 <- (1-u)^(1/v)
  u <- u1/(u1+u2)
  x <- logit(mu) + sqrt(2)*sigma*erfinv(2*u-1)
  q_ollltn <- 1/(1+exp(-x))
  return (q_ollltn)
}

sample <-function(mu, sigma, v) { # Returns a sample from OLLLTN using inverse transform
  u <- runif(1)
  return (Quantile(u,mu,sigma,v))
}


### Code to generate plot 2.1 (a) in the report
num_samples <- 10000 
mu <- c(0.15,0.25,0.50,0.60,0.80)
sigma <- c(0.25,0.30,0.50,0.20,0.20)
v <- c(0.20,0.20,0.20,0.20,0.25)
x1.vec <- rep(NA, num_samples-1)
x2.vec <- rep(NA, num_samples-1)
x3.vec <- rep(NA, num_samples-1)
x4.vec <- rep(NA, num_samples-1)
x5.vec <- rep(NA, num_samples-1)
y1.vec <- rep(NA, num_samples-1)
y2.vec <- rep(NA, num_samples-1)
y3.vec <- rep(NA, num_samples-1)
y4.vec <- rep(NA, num_samples-1)
y5.vec <- rep(NA, num_samples-1)
strt <- 0
ed <- c(0.6,0.8,0.97,0.83,0.98)
diff1 <- (ed[1]-strt)/num_samples
diff2 <- (ed[2]-strt)/num_samples
diff3 <- (ed[3]-strt)/num_samples
diff4 <- (ed[4]-strt)/num_samples
diff5 <- (ed[5]-strt)/num_samples
for(i in 2:num_samples){ 
  x1.vec[i-1] <- strt + (i-1)*diff1
  x2.vec[i-1] <- strt + (i-1)*diff2
  x3.vec[i-1] <- strt + (i-1)*diff3
  x4.vec[i-1] <- strt + (i-1)*diff4
  x5.vec[i-1] <- strt + (i-1)*diff5
} 

for(i in 2:num_samples){
  y1.vec[i-1] <- Pdf(mu[1], sigma[1], v[1], x1.vec[i-1])
  y2.vec[i-1] <- Pdf(mu[2], sigma[2], v[2], x2.vec[i-1])
  y3.vec[i-1] <- Pdf(mu[3], sigma[3], v[3], x3.vec[i-1])
  y4.vec[i-1] <- Pdf(mu[4], sigma[4], v[4], x4.vec[i-1])
  y5.vec[i-1] <- Pdf(mu[5], sigma[5], v[5], x5.vec[i-1])
}
pdf("Density1.pdf") #Graph is saved to this location
plot(x1.vec, y1.vec, type="l", col="navy blue", ylab="f(y)",xlab="y",lty=1, lwd=3.0, xlim=c(0,1.0), ylim=c(0,8))
points(x2.vec, y2.vec, col="dark red", type="l",lwd=3.0)
lines(x2.vec, y2.vec, col="dark red",lty=1, lwd=3.0,type="l")
points(x3.vec, y3.vec, col="orange", type="l",lwd=3.0)
lines(x3.vec, y3.vec, col="orange", lty=1,lwd=3.0,type="l")
points(x4.vec, y4.vec, col="green", type="l",lwd=3.0)
lines(x4.vec, y4.vec, col="green", lty=1,lwd=3.0,type="l")
points(x5.vec, y5.vec, col="dark green", type="l",lwd=3.0)
lines(x5.vec, y5.vec, col="dark green", lty=1,lwd=3.0,type="l")
legend(0.5,8,legend=c("mu=0.15,sigma=0.25,v=0.20","mu=0.25,sigma=0.30,v=0.20","mu=0.50,sigma=0.50,v=0.20","mu=0.60,sigma=0.20,v=0.20","mu=0.80,sigma=0.20,v=0.25"), col=c("navy blue","dark red","orange","green","dark green"),lty=c(1,1,1,1,1), ncol=1)
dev.off()

##########

### Plot to generate figure 2.1 (b) from report
num_samples <- 10000
mu <- c(0.5,0.5,0.5,0.5,0.5)
sigma <- c(0.5,0.5,0.5,0.5,0.5)
v <- c(0.50,0.45,0.35,0.25,0.15)
x1.vec <- rep(NA, num_samples-1)
x2.vec <- rep(NA, num_samples-1)
x3.vec <- rep(NA, num_samples-1)
x4.vec <- rep(NA, num_samples-1)
x5.vec <- rep(NA, num_samples-1)
y1.vec <- rep(NA, num_samples-1)
y2.vec <- rep(NA, num_samples-1)
y3.vec <- rep(NA, num_samples-1)
y4.vec <- rep(NA, num_samples-1)
y5.vec <- rep(NA, num_samples-1)
strt <- 0
ed <- c(1.0,1.0,1.0,1.0,0.983)
diff1 <- (ed[1]-strt)/num_samples
diff2 <- (ed[2]-strt)/num_samples
diff3 <- (ed[3]-strt)/num_samples
diff4 <- (ed[4]-strt)/num_samples
diff5 <- (ed[5]-strt)/num_samples
for(i in 2:num_samples){ 
  x1.vec[i-1] <- strt + (i-1)*diff1
  x2.vec[i-1] <- strt + (i-1)*diff2
  x3.vec[i-1] <- strt + (i-1)*diff3
  x4.vec[i-1] <- strt + (i-1)*diff4
  x5.vec[i-1] <- strt + (i-1)*diff5
} 

for(i in 2:num_samples){
  y1.vec[i-1] <- Pdf(mu[1], sigma[1], v[1], x1.vec[i-1])
  y2.vec[i-1] <- Pdf(mu[2], sigma[2], v[2], x2.vec[i-1])
  y3.vec[i-1] <- Pdf(mu[3], sigma[3], v[3], x3.vec[i-1])
  y4.vec[i-1] <- Pdf(mu[4], sigma[4], v[4], x4.vec[i-1])
  y5.vec[i-1] <- Pdf(mu[5], sigma[5], v[5], x5.vec[i-1])
}
pdf("Density2.pdf") #plot is saved to this location
plot(x1.vec, y1.vec, type="l", col="navy blue", ylab="f(y)",xlab="y",lty=1, lwd=3.0, ylim=c(0,3))
points(x2.vec, y2.vec, col="dark red", type="l",lwd=3.0)
lines(x2.vec, y2.vec, col="dark red",lty=1, lwd=3.0,type="l")
points(x3.vec, y3.vec, col="orange", type="l",lwd=3.0)
lines(x3.vec, y3.vec, col="orange", lty=1,lwd=3.0,type="l")
points(x4.vec, y4.vec, col="green", type="l",lwd=3.0)
lines(x4.vec, y4.vec, col="green", lty=1,lwd=3.0,type="l")
points(x5.vec, y5.vec, col="dark green", type="l",lwd=3.0)
lines(x5.vec, y5.vec, col="dark green", lty=1,lwd=3.0,type="l")
legend(0.5,3,legend=c("mu=0.50,sigma=0.50,v=0.50","mu=0.50,sigma=0.50,v=0.45","mu=0.50,sigma=0.50,v=0.35","mu=0.50,sigma=0.50,v=0.25","mu=0.50,sigma=0.50,v=0.15"), col=c("navy blue","dark red","orange","green","dark green"),lty=c(1,1,1,1,1), ncol=1)
dev.off()

#################################################

################## Section 3: Properties of OLLLTN ####################

skewness <-function(mu, sigma, v) { # Returns skewness for OLLLTN(mu,sigma,nu)
  num<- Quantile(1/4,mu,sigma,v) + Quantile(3/4,mu,sigma,v) - 2*Quantile(1/2,mu,sigma,v)
  den <- Quantile(3/4, mu, sigma, v) - Quantile(1/4, mu, sigma, v)
  return (num/den)
}
skewness_nu <-function(mu, sigma) { # Returns skewness for OLLLTN(mu,sigma,0.1)
  v=0.1
  num<- Quantile(1/4,mu,sigma,v) + Quantile(3/4,mu,sigma,v) - 2*Quantile(1/2,mu,sigma,v)
  den <- Quantile(3/4, mu, sigma, v) - Quantile(1/4, mu, sigma, v)
  return (num/den)
}
skewness_mu <-function(sigma, v) { # Returns skewness for OLLLTN(0.1,sigma,nu)
  mu=0.1
  num<- Quantile(1/4,mu,sigma,v) + Quantile(3/4,mu,sigma,v) - 2*Quantile(1/2,mu,sigma,v)
  den <- Quantile(3/4, mu, sigma, v) - Quantile(1/4, mu, sigma, v)
  return (num/den)
}
skewness_sigma <-function(mu, v) { # Returns skewness for OLLLTN(mu,0.1,nu)
  sigma=0.1
  num<- Quantile(1/4,mu,sigma,v) + Quantile(3/4,mu,sigma,v) - 2*Quantile(1/2,mu,sigma,v)
  den <- Quantile(3/4, mu, sigma, v) - Quantile(1/4, mu, sigma, v)
  return (num/den)
}
kurtosis <- function(mu, sigma, v) { # Returns kurtosis for OLLLTN(mu,sigma,nu)
  num <- Quantile(7/8, mu, sigma , v) - Quantile(5/8, mu, sigma , v) +Quantile(3/8, mu, sigma , v) -Quantile(1/8, mu, sigma , v)
  den <- Quantile(6/8, mu, sigma , v) - Quantile(2/8, mu, sigma , v)
  return (num/den)
}
kurtosis_nu <- function(mu, sigma) { # Returns kurtosis for OLLLTN(mu,sigma,0.1)
  v = 0.1
  num <- Quantile(7/8, mu, sigma , v) - Quantile(5/8, mu, sigma , v) +Quantile(3/8, mu, sigma , v) -Quantile(1/8, mu, sigma , v)
  den <- Quantile(6/8, mu, sigma , v) - Quantile(2/8, mu, sigma , v)
  return (num/den)
}
kurtosis_mu <- function(sigma, v) { # Returns kurtosis for OLLLTN(0.1,sigma,nu)
  mu = 0.1
  num <- Quantile(7/8, mu, sigma , v) - Quantile(5/8, mu, sigma , v) +Quantile(3/8, mu, sigma , v) -Quantile(1/8, mu, sigma , v)
  den <- Quantile(6/8, mu, sigma , v) - Quantile(2/8, mu, sigma , v)
  return (num/den)
}
kurtosis_sigma <- function(mu, v) { # Returns kurtosis for OLLLTN(mu,0.1,nu)
  sigma = 0.1
  num <- Quantile(7/8, mu, sigma , v) - Quantile(5/8, mu, sigma , v) +Quantile(3/8, mu, sigma , v) -Quantile(1/8, mu, sigma , v)
  den <- Quantile(6/8, mu, sigma , v) - Quantile(2/8, mu, sigma , v)
  return (num/den)
}
moment <-function(k, mu, sigma ,v ) { # Returns the kth moment for OLLLTN(mu,sigma,nu)
  #Returns the kth moment
  integrand <- function(u) {
    return (Quantile(u,mu,sigma,v)^k)
  }
  return (integrate(integrand, 0, 1)$value)
}


### Code to generate data for plots in section 3.2

num_val<-1000
x_p <- y_p <- seq(0.1, 1, length= num_val+1)
z_kurt_1 <- outer(x_p, y_p, kurtosis_nu)
z_kurt_2 <- outer(x_p, y_p, kurtosis_mu)
z_kurt_3 <- outer(x_p, y_p, kurtosis_sigma)

z_skew_1 <- outer(x_p, y_p, skewness_nu)
z_skew_2 <- outer(x_p, y_p, skewness_mu)
z_skew_3 <- outer(x_p, y_p, skewness_sigma)

write.csv(z_kurt_1,"z_kurt_nu.csv") 
write.csv(z_kurt_2,"z_kurt_mu.csv")
write.csv(z_kurt_3,"z_kurt_sigma.csv")

write.csv(z_skew_1,"z_skew_nu.csv")
write.csv(z_skew_2,"z_skew_mu.csv")
write.csv(z_skew_3,"z_skew_sigma.csv")


### Code to generate figures in section 3.3 of reports
num_samples <- 10000
mu <- c(0.30,0.50,0.80)
sigma <- c(0.30,0.50,0.25)
v <- c(0.30,0.25,0.21)
x1.vec <- rep(NA, num_samples-1)
x2.vec <- rep(NA, num_samples-1)
x3.vec <- rep(NA, num_samples-1)
y1.vec <- rep(NA, num_samples-1)
y2.vec <- rep(NA, num_samples-1)
y3.vec <- rep(NA, num_samples-1)
r1.vec <- rep(NA, num_samples-1)
r2.vec <- rep(NA, num_samples-1)
r3.vec <- rep(NA, num_samples-1)
strt <- 0
ed <- c(0.85,1.0,1.0)
diff1 <- (ed[1]-strt)/num_samples
diff2 <- (ed[2]-strt)/num_samples
diff3 <- (ed[3]-strt)/num_samples
for(i in 2:num_samples){ 
  x1.vec[i-1] <- strt + (i-1)*diff1
  x2.vec[i-1] <- strt + (i-1)*diff2
  x3.vec[i-1] <- strt + (i-1)*diff3
} 

for(i in 2:num_samples){
  y1.vec[i-1] <- Pdf(mu[1], sigma[1], v[1], x1.vec[i-1])*200
  y2.vec[i-1] <- Pdf(mu[2], sigma[2], v[2], x2.vec[i-1])*200
  y3.vec[i-1] <- Pdf(mu[3], sigma[3], v[3], x3.vec[i-1])*200
  r1.vec[i-1] <- sample(mu[1], sigma[1], v[1])
  r2.vec[i-1] <- sample(mu[2], sigma[2], v[2])
  r3.vec[i-1] <- sample(mu[3], sigma[3], v[3])
}
pdf("Density3a.pdf")
plot(x1.vec, y1.vec, type="l", col="navy blue", ylab="Density",xlab="y",lty=1, lwd=3.0, ylim=c(0,600))
hist(r1.vec,breaks=50,col=rgb(0.9,0.9,0.9,0.2),add=TRUE)
legend(0.4,600,legend=c("mu=0.30,sigma=0.30,v=0.30"), col=c("navy blue"),lty=c(1), ncol=1)
dev.off()
pdf("Density3b.pdf")
plot(x2.vec, y2.vec, type="l", col="dark red", ylab="Density",xlab="y",lty=1, lwd=3.0, ylim=c(0,400))
hist(r2.vec,breaks=50,col=rgb(0.9,0.9,0.9,0.2),add=TRUE)
legend(0.5,400,legend=c("mu=0.50,sigma=0.50,v=0.25"), col=c("dark red"),lty=c(1), ncol=1)
dev.off()
pdf("Density3c.pdf")
plot(x3.vec, y3.vec, type="l", col="dark green", ylab="Density",xlab="y",lty=1, lwd=3.0, ylim=c(0,1000))
hist(r3.vec,breaks=50,col=rgb(0.9,0.9,0.9,0.2),add=TRUE)
legend(0.0,1000,legend=c("mu=0.80,sigma=0.25,v=0.21"), col=c("dark green"),lty=c(1), ncol=1)
dev.off()


#################################################


#################### Section 4: OLLLTN Regression #######################

########################################################################


################## Section 5: Monte Carlo Simulation ####################

library(gamlss)
my_erf <- function(mu,sigma,y) { #erf function while accounting for precision errors

  x = (logit(y)-logit(mu))/(sqrt(2)*sigma)
  t = rep(0, length(y))
  for (i in 1:length(t)) {
    t[i] = erf(x[i])
  }
  t[t==1] = 0.99999999 # values that were rounded off to 1 should not be 1
  t[t==0] = 0.00000001 # values that were rounded off to 0 should not be 0
  t[t==-1] = -0.99999999 # values that were rounded off to -1 should not be -1
  return (t)
  
}

qOLLLTN <- function(p, mu = 0.1, sigma = 0.1, nu = 0.1, lower.tail=TRUE, log.p = FALSE){ 
  # quantile function for OLLLTN using gamlss (this was given in the paper)
  if (any(mu<0) | any(mu>1)) stop(paste("mu must be between 0 and 1", "nn", ""))
  if (any(sigma<=0)) {
    stop(paste("sigma must be positive","nn",""))
  }
  if (any(nu<=0)) {
    stop(paste("nu must be positive","nn",""))
  }
  if (any(p < 0)| any(p > 1)) stop(paste("p must be between 0 and 1", "nn", "")) 
  if (log.p==TRUE) p <- exp(p) else p <- p
  if (lower.tail==TRUE) p <- p else p <- 1-p
  u <- (p^(1/nu))/((1-p)^(1/nu)+p^(1/nu))
  q <- qLOGITNO(p = u, mu = mu, sigma = sigma) #returns sample from Logit Normal distribution using gamlss
  q
}

sample <-function(mu, sigma, v) { # returns sample from OLLLTN 
  u <- runif(1)
  return (qOLLLTN(u,mu,sigma,v))
}

ollltn_optim_fit <- function(n_epochs,mu,sigma,nu) { 
  # finds best fit parameters mu,sigma,nu for given y and returns the following:
  # 1. Average Bias over n_epochs
  # 2. Average Estimate over n_epochs
  # 3. Average Mean Square Error over n_epochs
  mle <- function(theta) { 
    # Returns the value of the negative log likelihood of OLLLTN for a given set of parameters
    mu <- theta[1]
    sigma <- theta[2]
    nu <- theta[3]
    mu <-exp(mu)/(1+exp(mu)) # link function
    sigma <- exp(sigma)
    nu <- exp(nu)
    l <- 0
    l <- l + n*log(nu)
    l<- l - n*log(sigma)
    l <- l - sum(((logit(y)-logit(mu))/sigma)^2)/2
    G <- (1/2)*(1+my_erf(mu,sigma,y))
    G[G==1] <- 1-1e-7 # Accounting for precision errors
    G[G==0] <- 1e-7
    l <- l+ (nu-1)*sum(log(G))
    l<- l+ (nu-1)*sum(log(1-G))
    l <- l - 2*sum(log(G^nu + (1-G)^nu))
    return(-l)
  }
  n <- 1000 # Number of Samples to generate
  bias <- c(0,0,0)
  mse <- c(0,0,0)
  j <- 0
  while(j<n_epochs) {
    print(j)
    y<-rep(0, n)
    for (i in 1:n) { 
      y[i] <- sample(mu,sigma,nu)
    }
    y[y==1] <- 1-1e-7 #Accounts for precision errors
    y[y==0] <- 1e-7
    skip_to_next <- FALSE
    tryCatch(fit <- optim(c(0.5, 0.5,0.5), mle, method = "L-BFGS-B",control=list(parscale=c(0.5,0.5,0.5))), error = function(e) { skip_to_next <<- TRUE})
    if(skip_to_next) { # Skip loop if optim does not converge
      next 
    }  
    j <- j+1
    mu_estimate <- fit$par[1]
    sigma_estimate <- fit$par[2]
    nu_estimate <- fit$par[3]
    mu_estimate <- exp(mu_estimate)/(1+exp(mu_estimate))
    sigma_estimate <- exp(sigma_estimate)
    nu_estimate <- exp(nu_estimate)
    bias <- bias + c(mu_estimate,sigma_estimate,nu_estimate)
    mse <- mse + c((mu_estimate-mu)^2,(sigma_estimate-sigma)^2,(nu_estimate-nu)^2)
  }
  
  avg_estimate <- bias/n_epochs
  bias <- c(mu,sigma,nu) - avg_estimate
  mse <- mse/n_epochs
  mse <- sqrt(mse)
  return (c(bias,avg_estimate, mse))
}
theta <- ollltn_optim_fit(10,0.8,0.25,0.21)
print(paste("Average Estimate of mu:",theta[4]))
print(paste("Average Estimate of sigma:",theta[5]))
print(paste("Average Estimate of nu:",theta[6]))
print(paste("Average Bias of mu:",theta[1]))
print(paste("Average Bias of sigma:",theta[2]))
print(paste("Average Bias of nu:",theta[3]))
print(paste("Average MSE of mu:",theta[7]))
print(paste("Average MSE of sigma:",theta[8]))
print(paste("Average MSE of nu:",theta[9]))



##############################################################################


################## Section 6: OLLLTN Regression Simulation ####################

ollltn_log_likelihood_regression <- function(theta) { 
  # Returns the negative log likelihood of OLLLTN for given values of Beta(ij) and nu
  b10 <- theta[1]
  b11 <- theta[2]
  b12 <- theta[3]
  b20 <- theta[4]
  b21 <- theta[5]
  b22 <- theta[6]
  nu <-theta[7]
  nu <- exp(nu) / (1+exp(nu)) #link function
  mu <- rep(0,n)
  sigma <- rep(0,n)
  mu <- exp(b10 + b11*x1+b12*x2) / (1+exp(b10 + b11*x1+b12*x2))
  sigma <- exp(b20 + b21*x1+b22*x2)

  G <- (1/2)*(1+my_erf(mu,sigma,y))
  G[G==1] <- 1-1e-7 # to account for precision
  G[G==0] <- 1e-7
  l<- 0
  l <- l + n*log(nu)
  l<- l - sum(log(sigma))
  l <- l - sum(((logit(y)-logit(mu))/sigma)^2)/2
  l <- l+ (nu-1)*sum(log(G))
  l<- l+ (nu-1)*sum(log(1-G))
  l <- l - 2*sum(log(G^nu + (1-G)^nu))

  return(-l)
}
logit_log_likelihood <- function(theta) { 
  # Returns the negative log likelihood of logit normal distribution for given parameters
  b10 <- theta[1]
  b11 <- theta[2]
  b12 <- theta[3]
  b20 <- theta[4]
  b21 <- theta[5]
  b22 <- theta[6]

  mu <- rep(0,n)
  sigma <- rep(0,n)
  mu <- exp(b10 + b11*x1+b12*x2) / (1+exp(b10 + b11*x1+b12*x2))
  sigma <- exp(b20 + b21*x1+b22*x2)
  
  l<- 0
  l<- l - sum(log(sigma))
  l <- l - sum(((logit(y)-logit(mu))/sigma)^2)/2
  return(-l)
}
n<-5000 # Number of samples
# Set initial value of Beta(ij) and nu
b10 <- 0.7 
b11 <- 0.1
b12 <- 0.9
b20 <- 0.3
b21 <- 0.5
b22 <- 0.2
nu <-0.3
n_epochs <- 5
theta <- c(0.7,0.1,0.9,0.3,0.5,0.2,0.3)

#Generate synthetic data
x1 <- rnorm(n,0,1)
x2 <- rbinom(n,1,0.5)


mu <- rep(0,n)
sigma <- rep(0,n)
y <- rep(0,n)
bias <- rep(0,7)
mse <- rep(0,7)
j <- 0
while(j<n_epochs) { #Carry out regression n_epochs time
  print(j)
  for (i in 1:n) { # Generate samples from OLLLTN(mu[i],sigma[i],nu[i])
    mu[i] <- exp(b10+b11*x1[i]+b12*x2[i]) / (1+exp(b10+b11*x1[i]+b12*x2[i]))
    sigma[i] <- exp(b20 + b21*x1[i]+b22*x2[i])
    y[i] <- sample(mu[i],sigma[i],nu)
  }
  y[abs(y-1)<(0.0000000000001)] <- 0.999999999 #for precision errors
  y[abs(y)<(0.0000000000001)] <- 0.000000001
  starting_vals_logit = rep(0.5,6) #Starting values
  logit_pars<-optim(starting_vals_logit, logit_log_likelihood, method = "L-BFGS-B",, control=list(parscale=starting_vals_logit))$par
  starting_vals <- c(logit_pars,1) #Use parameters from logit regression as starting values
  lower = c(rep(0.05,7))
  upper = c(rep(0.95,6),1.001)
  params <- optim(starting_vals, ollltn_log_likelihood_regression, method = "L-BFGS-B", control=list(parscale=starting_vals))$par
  j <- j+1
  params[7] <- exp(params[7])/(1+exp(params[7]))
  print(params)
  bias <- bias + params
  mse <- mse + (params-theta)^2
}
avg_estimate <- bias/n_epochs #AE for all parameters
bias <- theta - avg_estimate #Average Bias for all parameters
mse <- mse/n_epochs
mse <- sqrt(mse) #Average MSE for all parameters

#################################################


################## Section 7: Applications ####################

# Here we fit the data using the gamlss package
# First we define the OLLLTN family gamlss object using the already existing Logit Normal Family, then fit the data
oddLTN <- OLLLTN <- function (mu.link = "logit", sigma.link = "log", nu.link = "log") 
{
  mstats <- checklink("mu.link", "oddLTN", substitute(mu.link), 
                      c("logit", "probit", "cloglog", "cauchit", "log", "own"))
  dstats <- checklink("sigma.link", "oddLTN", substitute(sigma.link), 
                      c("inverse", "log", "identity"))
  vstats <- checklink("nu.link", "oddLTN", substitute(nu.link), 
                      c("log", "own","identity"))
  structure(list(family = c("OLLLTN", "oddLTN", "identity"), 
                 parameters = list(mu = TRUE, sigma = TRUE, nu = TRUE), 
                 nopar = 3, 
                 type = "Continuous",
                 mu.link = as.character(substitute(mu.link)), 
                 sigma.link = as.character(substitute(sigma.link)),
                 nu.link = as.character(substitute(nu.link)), 
                 mu.linkfun = mstats$linkfun, 
                 sigma.linkfun = dstats$linkfun, 
                 nu.linkfun = vstats$linkfun, 
                 mu.linkinv = mstats$linkinv, 
                 sigma.linkinv = dstats$linkinv, 
                 nu.linkinv = vstats$linkinv,
                 mu.dr = mstats$mu.eta, 
                 sigma.dr = dstats$mu.eta, 
                 nu.dr = vstats$mu.eta, 
                 dldm = function(y, mu, sigma, nu) {
                   e1 <- expression(-((log(y/(1-y))-log(mu/(1-mu)))^2)/(2*sigma^2) + log(nu))
                   e2<- expression((nu-1)*log(pnorm((log(y/(1-y))-log(mu/(1-mu)))/sigma)) + (nu-1)*log(pnorm((log((1-y)/y)-log(mu/(1-mu)))/sigma)))
                   e3 <- expression(-2*log((pnorm((log(y/(1-y))-log(mu/(1-mu)))/sigma))^nu+(pnorm((log((1-y)/y)-log(mu/(1-mu)))/sigma))^nu))
                   dldm <- eval(D(e1,'mu')) + eval(D(e2, 'mu')) + eval(D(e3,'mu'))
                   dldm
                 }, 
                 d2ldm2 = function(y, mu, sigma, nu) {
                   e1 <- expression(-((log(y/(1-y))-log(mu/(1-mu)))^2)/(2*sigma^2) + log(nu))
                   e2<- expression((nu-1)*log(pnorm((log(y/(1-y))-log(mu/(1-mu)))/sigma)) + (nu-1)*log(pnorm((log((1-y)/y)-log(mu/(1-mu)))/sigma)))
                   e3 <- expression(-2*log((pnorm((log(y/(1-y))-log(mu/(1-mu)))/sigma))^nu+(pnorm((log((1-y)/y)-log(mu/(1-mu)))/sigma))^nu))
                   e1 <- D(e1, 'mu')
                   e2 <- D(e2,'mu')
                   e3 <- D(e3, 'mu')
                   d2ldm2 <- eval(D(e1,'mu')) + eval(D(e2, 'mu')) + eval(D(e3,'mu'))
                   d2ldm2
                 },
                 dldd = function(y, mu, sigma, nu) {
                   e1 <- expression(-((log(y/(1-y))-log(mu/(1-mu)))^2)/(2*sigma^2) + log(nu))
                   e2<- expression((nu-1)*log(pnorm((log(y/(1-y))-log(mu/(1-mu)))/sigma)) + (nu-1)*log(pnorm((log((1-y)/y)-log(mu/(1-mu)))/sigma)))
                   e3 <- expression(-2*log((pnorm((log(y/(1-y))-log(mu/(1-mu)))/sigma))^nu+(pnorm((log((1-y)/y)-log(mu/(1-mu)))/sigma))^nu))
                   d2ldm2 <- eval(D(e1,'sigma')) + eval(D(e2, 'sigma')) + eval(D(e3,'sigma'))
                   d2ldm2
                 },
                 d2ldd2 = function(y, mu, sigma, nu) {
                   e1 <- expression(-((log(y/(1-y))-log(mu/(1-mu)))^2)/(2*sigma^2) + log(nu))
                   e2<- expression((nu-1)*log(pnorm((log(y/(1-y))-log(mu/(1-mu)))/sigma)) + (nu-1)*log(pnorm((log((1-y)/y)-log(mu/(1-mu)))/sigma)))
                   e3 <- expression(-2*log((pnorm((log(y/(1-y))-log(mu/(1-mu)))/sigma))^nu+(pnorm((log((1-y)/y)-log(mu/(1-mu)))/sigma))^nu))
                   e1 <- D(e1, 'sigma')
                   e2 <- D(e2,'sigma')
                   e3 <- D(e3, 'sigma')
                   d2ldm2 <- eval(D(e1,'sigma')) + eval(D(e2, 'sigma')) + eval(D(e3,'sigma'))
                   d2ldm2
                   
                 }, 
                 dldv = function(y, mu, sigma, nu) {
                   e1 <- expression(-((log(y/(1-y))-log(mu/(1-mu)))^2)/(2*sigma^2) + log(nu))
                   e2<- expression((nu-1)*log(pnorm((log(y/(1-y))-log(mu/(1-mu)))/sigma)) + (nu-1)*log(pnorm((log((1-y)/y)-log(mu/(1-mu)))/sigma)))
                   e3 <- expression(-2*log((pnorm((log(y/(1-y))-log(mu/(1-mu)))/sigma))^nu+(pnorm((log((1-y)/y)-log(mu/(1-mu)))/sigma))^nu))
                   dldv <- eval(D(e1,'nu')) + eval(D(e2, 'nu')) + eval(D(e3,'nu'))
                   dldv
                 }, 
                 d2ldv2 = function(y, mu, sigma, nu) {
                   e1 <- expression(-((log(y/(1-y))-log(mu/(1-mu)))^2)/(2*sigma^2) + log(nu))
                   e2<- expression((nu-1)*log(pnorm((log(y/(1-y))-log(mu/(1-mu)))/sigma)) + (nu-1)*log(pnorm((log((1-y)/y)-log(mu/(1-mu)))/sigma)))
                   e3 <- expression(-2*log((pnorm((log(y/(1-y))-log(mu/(1-mu)))/sigma))^nu+(pnorm((log((1-y)/y)-log(mu/(1-mu)))/sigma))^nu))
                   e1 <- D(e1, 'nu')
                   e2 <- D(e2,'nu')
                   e3 <- D(e3, 'nu')
                   d2ldv2 <- eval(D(e1,'nu')) + eval(D(e2, 'nu')) + eval(D(e3,'nu'))
                   d2ldv2
                   
                 },
                 d2ldmdd = function(y, mu, sigma, nu) {
                   e1 <- expression(-((log(y/(1-y))-log(mu/(1-mu)))^2)/(2*sigma^2) + log(nu))
                   e2<- expression((nu-1)*log(pnorm((log(y/(1-y))-log(mu/(1-mu)))/sigma)) + (nu-1)*log(pnorm((log((1-y)/y)-log(mu/(1-mu)))/sigma)))
                   e3 <- expression(-2*log((pnorm((log(y/(1-y))-log(mu/(1-mu)))/sigma))^nu+(pnorm((log((1-y)/y)-log(mu/(1-mu)))/sigma))^nu))
                   e1 <- D(e1, 'sigma')
                   e2 <- D(e2,'sigma')
                   e3 <- D(e3, 'sigma')
                   d2ldmdd <- eval(D(e1,'mu')) + eval(D(e2, 'mu')) + eval(D(e3,'mu'))
                   d2ldmdd
                 }, 
                 d2ldmdv = function(y, mu, sigma, nu) {
                   e1 <- expression(-((log(y/(1-y))-log(mu/(1-mu)))^2)/(2*sigma^2) + log(nu))
                   e2<- expression((nu-1)*log(pnorm((log(y/(1-y))-log(mu/(1-mu)))/sigma)) + (nu-1)*log(pnorm((log((1-y)/y)-log(mu/(1-mu)))/sigma)))
                   e3 <- expression(-2*log((pnorm((log(y/(1-y))-log(mu/(1-mu)))/sigma))^nu+(pnorm((log((1-y)/y)-log(mu/(1-mu)))/sigma))^nu))
                   e1 <- D(e1, 'mu')
                   e2 <- D(e2,'mu')
                   e3 <- D(e3, 'mu')
                   d2ldmdv <- eval(D(e1,'nu')) + eval(D(e2, 'nu')) + eval(D(e3,'nu'))
                   d2ldmdv
                 }, 
                 d2ldddv = function(y, mu, sigma, nu) {
                   e1 <- expression(-((log(y/(1-y))-log(mu/(1-mu)))^2)/(2*sigma^2) + log(nu))
                   e2<- expression((nu-1)*log(pnorm((log(y/(1-y))-log(mu/(1-mu)))/sigma)) + (nu-1)*log(pnorm((log((1-y)/y)-log(mu/(1-mu)))/sigma)))
                   e3 <- expression(-2*log((pnorm((log(y/(1-y))-log(mu/(1-mu)))/sigma))^nu+(pnorm((log((1-y)/y)-log(mu/(1-mu)))/sigma))^nu))
                   e1 <- D(e1, 'sigma')
                   e2 <- D(e2,'sigma')
                   e3 <- D(e3, 'sigma')
                   d2ldddv <- eval(D(e1,'nu')) + eval(D(e2, 'nu')) + eval(D(e3,'nu'))
                   d2ldddv
                 },
                 G.dev.incr = function(y, mu, sigma, nu, ...) -2 * dOLLLTN(y, mu, sigma,nu, log = TRUE),
                 rqres = expression(rqres(pfun = "pOLLLTN", type = "Continuous", y = y, mu = mu, sigma = sigma, nu = nu)), 
                 mu.initial = expression({mu <- rep(0.5, length(y))}), 
                 sigma.initial = expression({sigma <- rep(sd(log(y/(1 - y))), length(y))}), 
                 nu.initial = expression({nu <- rep(1, length(y))}), 
                 mu.valid = function(mu) all(mu > 0 & mu < 1), 
                 sigma.valid = function(sigma) all(sigma > 0), 
                 nu.valid = function(nu) all(nu > 0), 
                 y.valid = function(y) all(y > 0 & y < 1)),
            class = c("gamlss.family", "family"))
}
dOLLLTN <- function(y,mu, sigma, nu, log = FALSE){ # returns pdf of OLLLTN(mu,sigma,nu)
  #y must be less than 1?
  if (any(mu<0) | any(mu>1)) stop(paste("mu must be between 0 and 1", "nn", ""))
  if (any(sigma<=0)) {
    stop(paste("sigma must be positive","nn",""))
  }
  if (any(nu<=0)) {
    stop(paste("nu must be positive","nn",""))
  }
  g <- pnorm((log(y/(1-y))-log(mu/(1-mu)))/sigma)
  g[g==1] = 1-1e-7
  g[g==0] = 1e-7
  d <- nu*(g*(1-g))^(nu-1) / (y*(1-y)*sqrt(2*pi)*sigma*exp(((log(y/(1-y))-log(mu/(1-mu)))^2)/(2*sigma^2))*((g^nu + (1-g)^nu)^2))
  if(log==TRUE) {
    d <- log(d)
  }
  d
}

rOLLLTN <-function(n,mu, sigma, nu) {
  # Returns n samples generated from OLLLTN(mu,sigma,nu)
  if (any(mu<0) | any(mu>1)) stop(paste("mu must be between 0 and 1", "nn", ""))
  if (any(sigma<=0)) {
    stop(paste("sigma must be positive","nn",""))
  }
  if (any(nu<=0)) {
    stop(paste("nu must be positive","nn",""))
  }
  u <- runif(n,0,1)
  s <- qOLLLTN(u,mu,sigma,nu)
  s[s==1] = 1-1e-7
  s[s==0] = 1e-7
  s
}

pOLLLTN <- function (q,mu,sigma,nu) {
  # Returns cdf of OLLLTN(mu,sigma,nu) at q
  if (any(mu<0) | any(mu>1)) stop(paste("mu must be between 0 and 1", "nn", ""))
  if (any(sigma<=0)) {
    stop(paste("sigma must be positive","nn",""))
  }
  if (any(nu<=0)) {
    stop(paste("nu must be positive","nn",""))
  }
  g <- pnorm((log(q/(1-q)) - log(mu/(1-mu)))/sigma)
  p <- g^nu / (g^nu + (1-g)^nu)
  p
}

#Load the data
data <- read.csv("/brazilian_cities_data.csv")
y<-data$HDI

fit <- gamlssML(y,OLLLTN)
#Our estimates
mu_estimate <- fit$mu 
sigma_estimate <- fit$sigma
nu_estimate <-fit$nu
#fit quality parameters
gD <- fit$G.deviance
AIC <-fit$aic
SBC <- fit$sbc

#fit for other distributions using pre defined families in gamlss
fit_LTN <- gamlssML(y,LOGITNO)
fit_BE <- gamlssML(y,BE)
fit_SIMPLEX <- gamlssML(y,SIMPLEX)

#Generate plots in section 7
num_samples <- 10000
mu <- c(0.6785535,0.6790461,0.6842778,0.6790925)
sigma <- c(0.1309002,0.7926546,0.3652149,0.1663487)
v <- c(0.2367922)
x1.vec <- rep(NA, num_samples-1)
y1.vec <- rep(NA, num_samples-1)
y2.vec <- rep(NA, num_samples-1)
y3.vec <- rep(NA, num_samples-1)
y4.vec <- rep(NA, num_samples-1)
strt <- 0
ed <- 0.99
diff <- (ed-strt)/num_samples
for(i in 2:num_samples){ 
  x1.vec[i-1] <- strt + (i-1)*diff
  y1.vec[i-1] <- Pdf(mu[1], sigma[1], v[1], x1.vec[i-1])*9
  y2.vec[i-1] <- dSIMPLEX(x1.vec[i-1],mu[1], sigma[1])*0.7
} 
pdf("HDI_hist.pdf")
plot(x1.vec, y1.vec, type="l", col="red", ylab="Density",xlab="y",lty=1, lwd=3.0, ylim=c(0,60),xlim=c(0.45,0.9))
points(x1.vec, y2.vec, col="navy blue", type="l",lwd=3.0)
lines(x1.vec, y2.vec, col="navy blue",lty=1, lwd=3.0,type="l")
hist(y,breaks=24,col=rgb(0.9,0.9,0.9,0.2),add=TRUE)
legend(0.5,40,legend=c("mu=0.6785535,sigma=0.1309002,v=0.2367922"), col=c("red"),lty=c(1), ncol=1)
dev.off()

###############################################################