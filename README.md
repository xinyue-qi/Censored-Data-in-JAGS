# Censored-Data-in-JAGS

### Objective:
To establish an automatic approach to specifying the correct deviance function in ``JAGS``, we propose a simple and generic alternative modeling strategy
for the analysis of censored outcomes.

## The default approach for censored data modeling in JAGS
```
model{ # Model 1
        for (o in 1:O){ # O is the number of observed cases;
              Y[o] ˜ f(theta[o]) # f need to be specified for JAGS
        }
        for (j in 1:J){ # J is the number of censored observations;
              # Left censoring (R=0): lim[j,] = c(cut[j], inf);
              # Right censoring (R=2): lim[j,] = c(-inf, cut[j]);
              # Interval censoring (R=1): lim[j,] = c(cut1[j], cut2[j]);
              R[j] ˜ dinterval(Y[O+j], lim[j,])
              Y[O+j] ˜ f(theta[O+j])
        }
        # prior for theta’s
}
```
## The proposed approach in JAGS

```
model{ # Model 2
	# block 1: fully-observed
	for (o in 1:O){
		Y[o] ˜ f(theta[o]) # f need to be specified for JAGS
	}
	# block 2: left/right censoring
	for (c in 1:C){
		Z1[c] ˜ dbern(p[c])
		p[c] <- F(cut[c], theta[O+c])
	}
	# block 3: interval censoring
	for (i in 1:I){
		Z2[i] ˜ dbern(p[i])
		p[C+i] <- F(cut2[i], theta[O+C+i]) - F(cut1[i], theta[O+C+i])
	}
	# prior for theta’s
}
```

## Illustrative Examples:

### 1. Survival data

``Data``
The survival data (`leukemia`) are openly available in the ``R`` package ``“survival”`` v3.2-11.

``JAGS Model``

```
# Default
model{
	for (o in 1:O){ 
		Y[o] ˜ dexp(lambda[o]) # observed
		lambda[o] <- exp(b0+b1*group[o])
	}

	for (j in 1:J){
		R[j] ˜ dinterval(Y[O+j],lim[j]) # right-censored
		Y[O+j] ˜ dexp(lambda[O+j])
		lambda[O+j] <- exp(b0+b1*group[O+j])
	}
	b0 ˜ dnorm(0, tau0) # tau0 fixed at 0.01
	b1 ˜ dnorm(0, tau1) # tau1 fixed at 0.01
}

# Proposed
model{
	for (o in 1:O){
		Y[o] ˜ dexp(lambda[o]) # observed
		lambda[o] <- exp(b0 + b1*group[o])
	}

	for (c in 1:C){
		Z[c] ˜ dbern(p[c]) # censoring status
		p[c] <- pexp(cut[c],lambda[O+c]) # cumulative exp. dist.
		lambda[O+c] <- exp(b0 + b1*group[O+c])
	}
	b0 ˜ dnorm(0, 0.01)
	b1 ˜ dnorm(0, 0.01)
}
```

### 2. Binomial data

``Data``
The binomial data were derived from the literature listed in doi:10.1001/jamaoncol.2019.0393 and available in this repository (`AEdata.csv`).

``JAGS Model``

```
# Model A
model{
	for (j in 1:J1){
       		Y[j] ~ dbin(theta[j], N[j])
	}
 
       	for (j in 1:J2){
       		Z[j] ~ dbern(p[j])
		p[j] <- pbin(cut[j], theta, N[j+J1])  #Y<=cut
	}

	theta ~ dbeta(1,1)
}

# Model B
model{
	for (j in 1:J1){
       		Y[j] ~ dbin(theta[j], N[j])
		theta[j] <- theta.group[group2[j]]
	} 
       	for (j in 1:J2){
       		Z[j] ~ dbern(p[j])
		p[j] <- pbin(cut[j], theta[j+J1], N[j+J1])  #Y<=cut
		theta[j+J1] <- theta.group[group2[j+J1]]
	}
	for (i in 1:n.group2){
		theta.group[i] ~ dbeta(1, 1)
	}
}

# Model C
model{
	for (j in 1:J1){
       		Y[j] ~ dbin(theta[j], N[j])
		theta[j] <- theta.drug[drug[j]]
	} 
       	for (j in 1:J2){
       		Z[j] ~ dbern(p[j])
		p[j] <- pbin(cut[j], theta[j+J1], N[j+J1])  #Y<=cut
		theta[j+J1] <- theta.drug[drug[j+J1]]
	}
	for (i in 1:n.drug){
		theta.drug[i] ~ dbeta(1, 1)
	}
}

# Model D
model{
	for (j in 1:J1){
       		Y[j] ~ dbin(theta[j], N[j])
		logit(theta[j]) <- theta.drug[drug[j]]
	} 
       	for (j in 1:J2){
       		Z[j] ~ dbern(p[j])
		p[j] <- pbin(cut[j], theta[j+J1], N[j+J1])  #Y<=cut
		logit(theta[j+J1]) <- theta.drug[drug[j+J1]]
	}
	for (i in 1:n.drug){
		theta.drug[i] <- mu+sigma.drug*sn.drug[i]
		sn.drug[i] ~ dnorm(0, 1)
	}
	mu ~ dt(0, 10, 1)
	sigma.drug ~ dt(0, a, 1)T(0,)
}

# Model E
model{
	for (j in 1:J1){
       		Y[j] ~ dbin(theta[j], N[j])
		cloglog(theta[j]) <- theta.drug[drug[j]]
	} 
       	for (j in 1:J2){
       		Z[j] ~ dbern(p[j])
		p[j] <- pbin(cut[j], theta[j+J1], N[j+J1])  #Y<=cut
		cloglog(theta[j+J1]) <- theta.drug[drug[j+J1]]
	}
	for (i in 1:n.drug){
		theta.drug[i] <- mu+sigma.drug*sn.drug[i]
		sn.drug[i] ~ dnorm(0, 1)
	}
	sigma.drug ~ dt(0, a, 1)T(0,) # A = 2.5 
	mu ~ dt(0, 10, 1)
}

# Model F
model{
	for (j in 1:J1){
       		Y[j] ~ dbin(theta[j], N[j])
		probit(theta[j]) <- theta.drug[drug[j]]
	} 
       	for (j in 1:J2){
       		Z[j] ~ dbern(p[j])
		p[j] <- pbin(cut[j], theta[j+J1], N[j+J1])  #Y<=cut
		probit(theta[j+J1]) <- theta.drug[drug[j+J1]]
	}
	for (i in 1:n.drug){
		theta.drug[i] <- mu+sigma.drug*sn.drug[i]
		sn.drug[i] ~ dnorm(0, 1)
	}
	mu ~ dt(0, 10, 1)
	sigma.drug ~ dt(0, a, 1)T(0,)
}

# Model G
model{
	for (j in 1:J1){
       		Y[j] ~ dbin(theta[j], N[j])
		theta[j] ~ dbeta(1,1)
	}
 
       	for (j in 1:J2){
       		Z[j] ~ dbern(p[j])
		p[j] <- pbin(cut[j], theta[j+J1], N[j+J1])  #Y<=cut
		theta[j+J1] ~ dbeta(1,1)
	}
}

```
