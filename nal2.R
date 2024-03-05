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
ind <- c(1)
s <- 1
for (i in 2:length(tabela$drzava)) {
  if (tabela$drzava[i]!=tabela$drzava[i-1]) { s <-s+1 }
  ind <- c(ind,s)
}
tabela$indeks <- ind
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
dump(c("samples","samples.centr"),"s1")
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
samplesList.centr <- runMCMC(Cmcmc.centr, niter = 12000, nburnin = 2000,
                             nchains = 3, inits = initsFunction)
chainsPlot(samplesList.centr, burnin=0)
dump(c("samplesList.centr"),"verige2")


tabela2 <- tabela %>%
  group_by(indeks) %>%
  summarise(povprecje=mean(pricakovanadoba),n=length(pricakovanadoba),
            varianca=var(pricakovanadoba))

m <- length(tabela2$indeks)
n <- tabela2$n
xMatrix <- matrix(NA, ncol=m,nrow=max(n))
for (j in 1:m) {
  xMatrix[1:n[j],j] <- tabela[tabela$indeks==j,]$solanje-mean(tabela[tabela$indeks==j,]$solanje)
}
yMatrix <- matrix(NA, ncol=m,nrow=max(n))
for (j in 1:m) {
  yMatrix[1:n[j],j] <- tabela[tabela$indeks==j,]$pricakovanadoba
}


code2 <- nimbleCode({
  mu ~ dnorm(0,sd=100)
  eta ~ dunif(0,100)
  sigma ~ dunif(0,100)
  beta ~ dnorm(0,sd=100)
  etaBeta ~ dunif(0,100)
  
  for (j in 1:m) {
    muGroups[j] ~ dnorm(mu, sd = eta)
    betaGroups[j] ~ dnorm(beta, sd = etaBeta)
    for (i in 1:n[j]) {
      y[i, j] ~ dnorm(muGroups[j] + betaGroups[j] * x[i, j], sd = sigma);
    }
  }
})

constants2 <- list(m = m, n = n)

inits2 <- list(mu=mean(tabela2$povprecje),
               eta = sd(tabela2$povprecje),
               sigma=mean(sqrt(tabela2$varianca)),
               muGroups=tabela2$povprecje,
               betaGroups=rep(0,m),
               beta=0,
               etaBeta=1)

data2 <- list(y=yMatrix,x=xMatrix)

Rmodel2 <- nimbleModel(code = code2, constants = constants2,
                       inits = inits2, data = data2) 
Rmodel2$initializeInfo()

conf2 <- configureMCMC(Rmodel2)
conf2$addMonitors('muGroups', 'betaGroups')   #dodamo shranjevanje mu_j in beta_j

Rmcmc2 <- buildMCMC(conf2)
Cmodel2 <- compileNimble(Rmodel2)
Cmcmc2 <- compileNimble(Rmcmc2, project = Cmodel2)
samples2 <- runMCMC(Cmcmc2, niter = 12000, nburnin = 2000)

dump(c("samples2"),"hier")

samplesSummary(samples2)[c(1, 103, 2, 104, 102, 105, 205),]
samplesPlot(samples2, var = c("beta","betaGroups[1]"))
samplesPlot(samples2, var = c("mu","muGroups[1]"))
samplesPlot(samples2, var = c("etaBeta","eta","sigma"))

initsFunction2 <- function(){
  list(mu = rnorm(1, mean = mean(tabela2$povprecje), sd = 10),
       eta = runif(1, min = 0, max = 10),
       sigma = runif(1, min = 0, max = 10),
       muGroups = rnorm(m, mean = tabela2$povprecje, sd = 10),
       betaGroups = rnorm(m),                 
       beta = rnorm(1),                       
       etaBeta = runif(1, min = 0, max = 10)) 
}

samplesList2 <- runMCMC(Cmcmc2, niter = 12000, nburnin = 2000,
                        nchains = 3, inits = initsFunction2)
dump(c("samplesList2"),"hierverige")
