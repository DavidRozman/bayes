library(readr)
library(dplyr)
library(tidyr)
library(coda)
library(nimble)
library(basicMCMCplots)

tabela <- read_csv("Life Expectancy Data.csv") %>%
  select("Country","Year", "Life expectancy", "Adult Mortality",
         "BMI", "GDP", "Schooling", "Status") %>%
  rename(drzava=Country, leto=Year, pricakovanadoba="Life expectancy",
         smrtnost="Adult Mortality", ITM=BMI,BDP=GDP,solanje=Schooling,
         status=Status) %>%
  drop_na()
tabela$status <- ifelse(tabela$status=="Developed",1,0)
View(tabela)

code <- nimbleCode({
  beta0 ~ dnorm(0, sd = 100)
  for(k in 1:p) {
    beta[k] ~ dnorm(0, sd = 100)
  }
  sigma ~ dunif(0, 100)
  for(i in 1:n) {
    y[i] ~ dnorm(beta0 + inprod(beta[1:p], x[i, 1:p]), sd = sigma)
  }
})

X <- subset(tabela, select = c("BDP","ITM","solanje","smrtnost","status"))
p <- ncol(X)

constants <- list(n = length(tabela$pricakovanadoba),
                  p = p,
                  x = X)

data <- list(y = tabela$pricakovanadoba)

inits <- list(beta0 = mean(tabela$pricakovanadoba),
              beta = rep(0, p),
              sigma = 1)

Rmodel <- nimbleModel(code, constants, data, inits)

conf <- configureMCMC(Rmodel)

Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)

Cmcmc <- compileNimble(Rmcmc, project=Cmodel)

samples <- runMCMC(Cmcmc, niter=12000,nburnin = 2000)
dump(c("samples"),"s1")
samplesSummary(samples)

samplesPlot(samples)
effectiveSize(samples)
cor(samples)

X.centr = apply(X, 2, function(y) y - mean(y))

constants.centr <- list(n = length(tabela$pricakovanadoba),
                        p = p,
                        x = X.centr)

Rmodel.centr <- nimbleModel(code, constants.centr, data, inits)
conf.centr <- configureMCMC(Rmodel.centr)
Rmcmc.centr <- buildMCMC(conf.centr)
Cmodel.centr <- compileNimble(Rmodel.centr)
Cmcmc.centr <- compileNimble(Rmcmc.centr, project = Cmodel.centr)
samples.centr <- runMCMC(Cmcmc.centr, niter = 12000, nburnin = 2000)

samplesPlot(samples)       #pred centriranjem
samplesPlot(samples.centr) #po centriranju

effectiveSize(samples.centr) #zelo dobro

cor(samples.centr) #veliko manjse kot prej

initsFunction <- function(){
  list(beta0 = rnorm(1),
       beta = rnorm(5),
       sigma = runif(1, min = 0, max = 10))
}

samplesList <- runMCMC(Cmcmc, niter = 12000, nburnin = 0,
                       nchains = 3, inits = initsFunction)
chainsPlot(samplesList, burnin=2000)
dump(c("samplesList"),"verige1")
