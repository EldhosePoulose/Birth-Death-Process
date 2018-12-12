## General Parameters
set.seed(123)

## birth and death process:
RUNS <- 10
T <- 500 # max time
alpha <- 0.1 # individual death rate
beta <- 0.1 # individual reproduction rate

for (run in seq(RUNS)){
  i = 1
  t = c(0)
  y = c(100)
  while(t[i] < T && y[i] > 0){
    i = i+1
    lambda = beta*y[i-1]
    mu = alpha*y[i-1]
    t[i] = t[i-1] + rexp(n=1, rate = mu+lambda) #sojourn time is exponentially distributed
    y[i] = y[i-1] + rbinom(n=1, prob = lambda/(lambda+ mu), size = 1) *2 -1
  }
  if(run == 1){
    plot(t, y, t='s', 
         col = run, 
         main = "Stochastic Birth and Death Process", 
         xlim = c(0, T+1), 
         ylim = c(0, 300),
         ylab = "Population Size",
         xlab = "Generations")
  } else {
    points(t, y, t='s', col = run)
  }
}



## sterile male insect control model
RUNS <- 50 # model repetitions
T <- 500 # max time
alpha <- 1 # individual death rate
beta <- 4 # individual reproduction rate
K <- 100 # carrying capacity
S <- 100 # number of sterile males

for (run in seq(RUNS)){
  i = 1
  t = c(0)
  y = c(20)
  while(t[i] < T && y[i] > 0){
    i = i+1
    if(y[i-1] < K){# Formula 5.15 from p.390
      lambda = beta * y[i-1] * y[i-1]/(y[i-1]+S) 
    } else {
      lambda = 0
    }
    mu = alpha * y[i-1]
    t[i] = t[i-1] + rexp(n=1, rate = mu+lambda) #sojourn time is exponentially distributed
    y[i] = y[i-1] + rbinom(n=1, prob = lambda/(lambda+ mu), size = 1) *2 -1
  }
  if(run == 1){
    plot(t, y, t='s', 
         col = run, 
         main = "Sterile Male Insect Control Model", 
         xlim = c(0, 10), 
         ylim = c(0, K),
         ylab = "Population Size",
         xlab = "Generations")
  } else {
    points(t, y, t='s', col = run)
  }
}


# Problem 5.2a
RUNS <- 9 # model repetitions
T <- 10 # max time
lambdas <- c(0,1,2,3,4,0)
mus <- c(0,4,3,2,1,0)

for (run in seq(RUNS)){
  i = 1
  t = c(0)
  y = c(2)
  while(t[i] < T && y[i] > 0 && y[i]<5){ # let T stay between 0 and 5
    i = i+1
    lambda <- lambdas[y[i-1]+1]
    mu <- mus[y[i-1]+1]
    t[i] = t[i-1] + rexp(n=1, rate = mu+lambda) #sojourn time is exponentially distributed
    y[i] = y[i-1] + rbinom(n=1, prob = lambda/(lambda+ mu), size = 1) *2 -1
  }
  if(run == 1){
    plot(t, y, t='l', 
         col = run, 
         main = "Sterile Male Insect Control Model, Prob 5.2a", 
         xlim = c(0, 2), 
         ylim = c(0, 5))
  } else {
    points(t, y, t='l', col = run)
  }
}

# Problem 5.2
lambda <- c(0,1,2,3,4,0)
mu <- c(0,4,3,2,1,0)
#transition probabilities:
#lambda=mu=0 for i=0 and 6 <- p[1]=p[6]=q[1]=q[6]=NA
#which is ok, because we don't need them
p <- lambda/(lambda+mu)
q <- 1-p
Q <- matrix(data = c(0, p[2], 0, 0,
                     q[3], 0, p[3], 0,
                     0, q[4], 0, p[4],
                     0, 0, q[5], 0),
            ncol = 4, nrow = 4, byrow = TRUE)
R <- matrix(data = c(q[2], 0,
                     0,0,
                     0,0,
                     0,p[5]),
            nrow = 4, ncol = 2, byrow = TRUE)
N = solve(diag(4)-Q) # fundamental Matrix
rates <- lambda + mu
absorbtion_time <- sum(N[2,]/rates[2:5])
B <- N%*%R

print(paste("absorbtion time", absorbtion_time))
print(paste("fundamental Matrix")); print(N)
print(paste("B")); print(B)
