# Survival Data

```
# load package
library(survival)

# alternative strategy
leuk.obs<-subset(leukemia,status==1)
leuk.cen<-subset(leukemia,status==0)
leuk.all<-rbind(leuk.obs,leuk.cen)
Y = leuk.obs$time
Z = leuk.cen$status     
cut = leuk.cen$time     # censoring time
J1 <-length(Y)
J2 <-length(Z)
group <- as.numeric(as.factor(leuk.all$x))

tau0 <- tau1 <- 0.01
data <-list(J1=J1, J2=J2, Y=Y, Z=Z, cut=cut, group=group, tau0=tau0, tau1=tau1)
inits <- function() {list(b0=rnorm(1),b1=rnorm(1))}
parameters <- c("lambda", "b0", "b1") 

library(rjags)
system.time(model.2 <- jags.model("modelfile.txt", data, inits, 
                                  n.chains=3, n.adapt=1000))
system.time(update(model.2,30000))
system.time(fit.2 <- coda.samples(model=model.2, variable.names=parameters,
                                  n.iter=30000, thin=3))

res2 <- round(summary(fit.2)[[2]][,c(3,1,5)],4)
```

# Binomial Data

```
# load study-level pneumonitis data
datal15 <- read.csv("AEdata.csv", header=T)   # 112 observations

# R is censoring variable 
# R=0 if left-censored (LC); R=1 if observed; R=2 if right-censored (RC)
# L is a cutoff if left-censored (note: L=-1 if fully observed)

data.obs <- datal15[!is.na(datal15$Y), ] # subset: 66 
data.mis <- datal15[is.na(datal15$Y), ]  # subset: 46 
data.all <- rbind(data.obs, data.mis)    # total: 112 

Z <- ifelse(data.mis$R==2, 0, 1)         # type of censoring (RC: Z=0, LC: Z=1)
cut <- data.mis$L                        # Left-censored: cut = L; 

Y <- as.numeric(data.obs$Y) # true number of AEs (observed)
Z <- as.numeric(Z)          # "missing" AEs
cut <- as.numeric(cut)
J1 <- length(Y)
J2 <- length(Z)
N <- as.numeric(data.all$No..of.treated.patients)  # number of pts

drug <- as.numeric(factor(data.all$Drug, 
                          level = c("Nivolumab", "Pembrolizumab",
                                    "Atezolizumab", "Avelumab",
                                    "Durvalumab"))) # 5 levels

n.drug <- length(table(drug))     # [1] 5

# Group = 1: PD1; Group = 2: PD-L1 
group2 <- as.numeric(as.factor(data.all$New.Drug.Category..target.))
n.group2 <- length(table(group2))

library(rjags)

# model A
data <- list (Y=Y, N=N, Z=Z, cut=cut, J1=J1, J2=J2)
inits <- function() {list (theta=rbeta(1,1,1))}
parameters <- c("theta")
# model B
data <- list (Y=Y, N=N, Z=Z, cut=cut, J1=J1, J2=J2, 
              group2=group2, n.group2=n.group2)
inits <- function() {list (theta.group=rep(0,2))}
parameters <- c("theta.group")
# model C
data <- list (Y=Y, N=N, Z=Z, cut=cut, J1=J1, J2=J2, 
              drug=drug, n.drug=n.drug)
inits <- function() {list (theta.drug=rep(0,5))}
parameters <- c("theta.drug")
# model D,E,F
data <- list (Y=Y, N=N, Z=Z, cut=cut, J1=J1, J2=J2, drug=drug, 
              n.drug=n.drug, a=a) # a = 0.16
inits <- function() {list (mu=rnorm(1))}
parameters <- c("mu.adj","p.drug.adj","sigma.drug")
# model G
data <- list (Y=Y, N=N, Z=Z, cut=cut, J1=J1, J2=J2)
inits <- function() {list (theta=rep(0,112))}
parameters <- c("theta")


model.3 <- jags.model ("modelfile.txt", data, inits,
                       n.chains=3, n.adapt=1000)
update(model.3,30000)
fit.3 <- coda.samples(model=model.3, variable.names=parameters,
                                  n.iter=30000, thin=3)
                                  
dic.pd <- dic.samples(model=model.3, n.iter=30000, type="pD"); dic.pd
dic.popt <- dic.samples(model=model.3, n.iter=30000, type="popt"); dic.popt

```
