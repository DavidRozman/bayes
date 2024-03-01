library(R2BayesX)

n <- 100
x1 <- runif(n,10,30)
x2 <- rnorm(n,80,10)
x3 <- rnorm(n,10,3)
x4 <- rnorm(n,5,2)
x5 <- runif(n,0,50)
y <- 20+4*x1+10*x2-0.2*x3-x4+rnorm(n,0,20)


plot(x4,y)


fit.lm <- lm(y~x1+x2+x3+x4+x5) 
summary(fit.lm)
b1 <- summary(fit.lm)$coefficients[2,1]
b2 <- summary(fit.lm)$coefficients[3,1]
b3 <- summary(fit.lm)$coefficients[4,1]
b4 <- summary(fit.lm)$coefficients[5,1]
napaka1 <- summary(fit.lm)$coefficients[2,2]
napaka2 <- summary(fit.lm)$coefficients[3,2]
napaka3 <- summary(fit.lm)$coefficients[4,2]
napaka4 <- summary(fit.lm)$coefficients[5,2]

fit.bayesx2 <- bayesx(y ~ x1+x2+x3+x4+x5, family = "gaussian", method = "MCMC")
summary(fit.bayesx)

bayes1 <- attr(fit.bayesx2$fixed.effects,"sample")[,2]
plot(bayes1,type="l")
hist(bayes1,prob=T)
lines(density(bayes1),col="red",lwd=2)

ba1<-fit.bayesx2$fixed.effects[2,1]
napaka1bayes <- summary(fit.lm)$coefficients[2,2]

NormalnaNapaka <- function(iteracije,varianca) {

b0lm <- c()    
b1lm <- c()
b2lm <- c()
b3lm <- c()
b4lm <- c()
b5lm <- c()
sn0 <- c()
sn1 <- c()
sn2 <- c()
sn3 <- c()
sn4 <- c()
sn5 <- c()
skn0 <- c()
skn1 <- c()
skn2 <- c()
skn3 <- c()
skn4 <- c()
skn5 <- c()

b0bayes <- c()
b1bayes <- c()
b2bayes <- c()
b3bayes <- c()
b4bayes <- c()
b5bayes <- c()
odklon0 <- c()
odklon1 <- c()
odklon2 <- c()
odklon3 <- c()
odklon4 <- c()
odklon5 <- c()
skn0b <- c()
skn1b <- c()
skn2b <- c()
skn3b <- c()
skn4b <- c()
skn5b <- c()


for (i in 1:iteracije) {
  n <- 100
  x1 <- runif(n,10,30)
  x2 <- rnorm(n,80,10)
  x3 <- rnorm(n,10,3)
  x4 <- rnorm(n,5,2)
  x5 <- runif(n,0,50)
  y <- 20+4*x1+10*x2-0.2*x3-x4+rnorm(n,0,varianca)
  
  fit.lm <- lm(y~x1+x2+x3+x4+x5)
  b0lm <- c(b0lm,summary(fit.lm)$coefficients[1,1])
  b1lm <- c(b1lm,summary(fit.lm)$coefficients[2,1])
  b2lm <- c(b2lm,summary(fit.lm)$coefficients[3,1])
  b3lm <- c(b3lm,summary(fit.lm)$coefficients[4,1])
  b4lm <- c(b4lm,summary(fit.lm)$coefficients[5,1])
  b5lm <- c(b5lm,summary(fit.lm)$coefficients[6,1])
  sn0 <- c(sn0,summary(fit.lm)$coefficients[1,2])
  sn1 <- c(sn1,summary(fit.lm)$coefficients[2,2])
  sn2 <- c(sn2,summary(fit.lm)$coefficients[3,2])
  sn3 <- c(sn3,summary(fit.lm)$coefficients[4,2])
  sn4 <- c(sn4,summary(fit.lm)$coefficients[5,2])
  sn5 <- c(sn5,summary(fit.lm)$coefficients[6,2])
  skn0 <- c(skn0,(20-summary(fit.lm)$coefficients[1,1])^2)
  skn1 <- c(skn1,(4-summary(fit.lm)$coefficients[2,1])^2)
  skn2 <- c(skn2,(10-summary(fit.lm)$coefficients[3,1])^2)
  skn3 <- c(skn3,(-0.2-summary(fit.lm)$coefficients[4,1])^2)
  skn4 <- c(skn4,(-1-summary(fit.lm)$coefficients[5,1])^2)
  skn5 <- c(skn5,(0-summary(fit.lm)$coefficients[6,1])^2)
  
  fit.bayesx <- bayesx(y ~ x1+x2+x3+x4+x5, family = "gaussian", method = "MCMC")
  b0bayes <- c(b0bayes,fit.bayesx$fixed.effects[1,1])
  b1bayes <- c(b1bayes,fit.bayesx$fixed.effects[2,1])
  b2bayes <- c(b2bayes,fit.bayesx$fixed.effects[3,1])
  b3bayes <- c(b3bayes,fit.bayesx$fixed.effects[4,1])
  b4bayes <- c(b4bayes,fit.bayesx$fixed.effects[5,1])
  b5bayes <- c(b5bayes,fit.bayesx$fixed.effects[6,1])
  odklon0 <- c(odklon0,fit.bayesx$fixed.effects[1,2])
  odklon1 <- c(odklon1,fit.bayesx$fixed.effects[2,2])
  odklon2 <- c(odklon2,fit.bayesx$fixed.effects[3,2])
  odklon3 <- c(odklon3,fit.bayesx$fixed.effects[4,2])
  odklon4 <- c(odklon4,fit.bayesx$fixed.effects[5,2])
  odklon5 <- c(odklon5,fit.bayesx$fixed.effects[6,2])
  skn0b <- c(skn0b,(20-fit.bayesx$fixed.effects[1,1])^2)
  skn1b <- c(skn1b,(4-fit.bayesx$fixed.effects[2,1])^2)
  skn2b <- c(skn2b,(10-fit.bayesx$fixed.effects[3,1])^2)
  skn3b <- c(skn3b,(-0.2-fit.bayesx$fixed.effects[4,1])^2)
  skn4b <- c(skn4b,(-1-fit.bayesx$fixed.effects[5,1])^2)
  skn5b <- c(skn5b,(0-fit.bayesx$fixed.effects[6,1])^2)
}
M <- matrix(0,nrow=6,ncol=6)
M[1,1] <- mean(b0lm)
M[2,1] <- mean(b1lm)
M[3,1] <- mean(b2lm)
M[4,1] <- mean(b3lm)
M[5,1] <- mean(b4lm)
M[6,1] <- mean(b5lm)

M[1,2] <- mean(sn0)
M[2,2] <- mean(sn1)
M[3,2] <- mean(sn2)
M[4,2] <- mean(sn3)
M[5,2] <- mean(sn4)
M[6,2] <- mean(sn5)

M[1,3] <- sqrt(mean(skn0))
M[2,3] <- sqrt(mean(skn1))
M[3,3] <- sqrt(mean(skn2))
M[4,3] <- sqrt(mean(skn3))
M[5,3] <- sqrt(mean(skn4))
M[6,3] <- sqrt(mean(skn5))

M[1,4] <- mean(b0bayes)
M[2,4] <- mean(b1bayes)
M[3,4] <- mean(b2bayes)
M[4,4] <- mean(b3bayes)
M[5,4] <- mean(b4bayes)
M[6,4] <- mean(b5bayes)

M[1,5] <-  mean(odklon0)
M[2,5] <-  mean(odklon1)
M[3,5] <-  mean(odklon2)
M[4,5] <-  mean(odklon3)
M[5,5] <-  mean(odklon4)
M[6,5] <-  mean(odklon5)

M[1,6] <- sqrt(mean(skn0b))
M[2,6] <- sqrt(mean(skn1b))
M[3,6] <- sqrt(mean(skn2b))
M[4,6] <- sqrt(mean(skn3b))
M[5,6] <- sqrt(mean(skn4b))
M[6,6] <- sqrt(mean(skn5b))

M
}


m0 <- NormalnaNapaka(100,50)
m1 <- NormalnaNapaka(100,50)
m2 <- NormalnaNapaka(100,50)
m3 <- NormalnaNapaka(100,50)
m4 <- NormalnaNapaka(100,50)
m5 <- NormalnaNapaka(100,50)
m6 <- NormalnaNapaka(100,50)
m7 <- NormalnaNapaka(100,50)
m8 <- NormalnaNapaka(100,50)
m9 <- NormalnaNapaka(100,50)

Varianca50 <- (m0+m1+m2+m3+m4+m5+m6+m7+m8+m9)/10

dump(c("Varianca50"),"sigma50")


n0 <- NormalnaNapaka(100,2)
n1 <- NormalnaNapaka(100,2)
n2 <- NormalnaNapaka(100,2)
n3 <- NormalnaNapaka(100,2)
n4 <- NormalnaNapaka(100,2)
n5 <- NormalnaNapaka(100,2)
n6 <- NormalnaNapaka(100,2)
n7 <- NormalnaNapaka(100,2)
n8 <- NormalnaNapaka(100,2)
n9 <- NormalnaNapaka(100,2)

Varianca2 <- (n0+n1+n2+n3+n4+n5+n6+n7+n8+n9)/10

dump(c("Varianca2"),"sigma2")
