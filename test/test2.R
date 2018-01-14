# Generate data
# parameters for data generating
<<<<<<< HEAD
p <- 10# number of total covariates
pt <- 5# number of covariates affecting treatment effects
py <- 4# number of covariates affecting outcomes but not treatment effects
=======
p <- 5# number of total covariates
pt <- 2# number of covariates affecting treatment effects
py <- 2# number of covariates affecting outcomes but not treatment effects
>>>>>>> iss01
asym <- .5 # whether treatment effects are distributed asymmetrically across treated and control
n <- 10000 # total size of the dataset
propens <- .5 #treatment probability
treatsize <- .5 # treatment effect size
levsize <- .3

# draw treatment W
w <- rbinom(n, 1, propens)
w_ <- 1-w

# draw covariates X
X <- matrix(rnorm(n * p), nrow = n, ncol = p)

# generate treatment effects as function of X
if (p<pt+py) print("error: p>=pt+py required")
tau <- 0
for (iii in 1:pt) {
  tau <- tau + treatsize*pmax(X[,iii],array(0,n))*(2*X[,iii])
}

# generate average value of outcomes
mu <- treatsize*rowSums(X[,1:pt]) + levsize*rowSums(X[,(pt+1):(pt+py)])

y <- vector(length = n)
y_ <- vector(length = n)
censor <- vector(length = n)
scale <- vector(length = n)
scale_ <- vector(length = n)

scale1 = 15
for (i in 1:n){
  scale[i] <- mu[i] + asym * w[i] * tau[i] + (asym-1) * (1-w[i]) * tau[i]
  scale_[i] <- mu[i] + asym * w_[i] * tau[i] + (asym-1) *(1-w_[i]) * tau[i]
  # y[i] <- runif(1,50,1 / exp(-scale1 - scale[i]))
  y[i] <- rweibull(1, shape = 6, scale = 1 / exp(-scale - scale[i]))
  y_[i] <- rweibull(1, shape = 6, scale = 1/ exp(-scale - scale_[i]))
  if (y[i] > 3650)
    y[i] = 3650
  if (y_[i] > 3650)
    y[i] = 3650
  censor[i] <- runif(1,365,3650)
}


# generate
time = pmin(y,censor)  #observed time is min of censored and true
event = as.numeric(time==y)   # set to 1 if event is observed
# event = rep(1,length(y))

f <- paste("x", 1:p, sep="", collapse="+")

name <- c( name,  "y", "w", "completeCase", "tau_true", "propensity")

tau_true <- (1-2*w)*(y_ - y)

ntr <- round(.5*n)

dataTrain <- data.frame(X[1:ntr,], y[1:ntr], w[1:ntr], event[1:ntr], tau_true[1:ntr], propens)

names(dataTrain)=name

cv.option.temp = "CT"


tree<- causalTree(as.formula(paste("y~",f)),
                    data=dataTrain, treatment=dataTrain$w,
                    split.Rule="survival", split.Honest=F, cv.option="CT", minsize = 100,
                    split.alpha = 1, cv.alpha = 1, xval=0, cp=0, propensity = dataTrain$propensity,
                    completeCase = dataTrain$completeCase)

tree3 <- causalTree(as.formula(paste("y~",f)),
                    data=dataTrain, treatment=dataTrain$w,
                    split.Rule="CT", split.Honest=F, cv.option="CT", minsize = 500,
                    split.alpha = 1, cv.alpha = 1, xval=0, cp=0, propensity = dataTrain$propensity,
                    completeCase = dataTrain$completeCase)

