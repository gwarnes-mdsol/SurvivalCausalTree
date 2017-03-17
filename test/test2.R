# Generate data
# parameters for data generating
p <- 10# number of total covariates
pt <- 4# number of covariates affecting treatment effects
py <- 6# number of covariates affecting outcomes but not treatment effects
asym <- .5 # whether treatment effects are distributed asymmetrically across treated and control
n <- 1000 # total size of the dataset
propens <- .5 #treatment probability
sd = 30
treatsize <- .5 # treatment effect size
treatsize <- .3
levsize <- 1
censorProb <- .3 # proportion of censoring observations

# draw treatment W
w <- rbinom(n, 1, propens)
w_ <- 1-w

# draw covariates X
X <- matrix(rnorm(n * p), nrow = n, ncol = p)

if (p<pt+py) print("error: p>=pt+py required")
tau <- 0
for (iii in 1:pt) {
  tau <- tau + treatsize*pmax(X[,iii],array(0,n))*(2*X[,iii])
}

# generate average value of outcomes
mu <- abs(treatsize*rowSums(X[,1:pt])+levsize*rowSums(X[,(pt+1):(pt+py)]))

# generate complete case indicator
lambdaT = .002 # baseline hazard
lambdaC = .300  # hazard of censoring

# generate outcomes as function of treatment status, mu, tau, and noise
y <- rweibull(n, shape=1, scale = lambdaT * exp( mu + asym * w * tau + (1-asym) * (1-w) * tau))
y_ <- rweibull(n, shape=1, scale = lambdaT * exp( mu + asym * w_ * tau + (1-asym) * (1-w_) * tau))

# generate
Censor = rweibull(n, shape=1, scale=lambdaC)   #censoring time
time = pmin(y,Censor)  #observed time is min of censored and true
event = as.numeric(time==y)   # set to 1 if event is observed

f <- ""
nextx <- ""
if (p>1) {
  for (ii in 1:(p-1)) {
    nextx <- paste("x",ii, sep="")
    if (ii==1) {name <- nextx}
    if (ii>1) {name <- c(name, nextx)}
    f <- paste(f, nextx, "+", sep="")
  }
  f <- paste(f, "x", ii+1, sep="")
} else if (p==1) {
  f <- "x1"
}

for (ii in 1:p) {
  nextx <- paste("x",ii, sep="")
  if (ii==1) {name <- nextx}
  if (ii>1) {name <- c(name, nextx)}
}

name <- c( name,  "y", "w", "completeCase", "tau_true")

tau_true <- (1-2*w)*(y_ - y)

ntr <- round(.333*n)
nest <- round(.333*n)
ntest <- n - ntr - nest

dataTrain <- data.frame(X[1:ntr,], y[1:ntr], w[1:ntr], event[1:ntr], tau_true[1:ntr])
dataEst <- data.frame(X[(ntr+1):(ntr+nest),], y[(ntr+1):(ntr+nest)], w[(ntr+1):(ntr+nest)], event[(ntr+1):(ntr+nest)], tau_true[(ntr+1):(ntr+nest)])
dataTest <- data.frame(X[(ntr+nest+1):n,], y[(ntr+nest+1):n], w[(ntr+nest+1):n], event[(ntr+nest+1):n], tau_true[(ntr+nest+1):n])

names(dataTrain)=name
names(dataEst)=name
names(dataTest)=name


cv.option.temp = "CT"


tree2 <- causalTree(as.formula(paste("y~",paste(f))),
                    data=dataTrain, treatment=dataTrain$w,
                    split.Rule="survival", split.Honest=F, cv.option="CT", minsize = 50,
                    split.alpha = 1, cv.alpha = 1, xval=0, cp=0, propensity = rep(.5,999), completeCase = dataTrain$completeCase)
