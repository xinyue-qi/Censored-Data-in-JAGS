# Censored-Data-in-JAGS

## The default approach
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
## The proposed approach

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
The survival data are openly available in the ``R`` package ``“survival”`` v3.2-11.

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
The binomial data were derived from the literature listed in doi:10.1001/jamaoncol.2019.0393 and available in the public domain.

``JAGS Model``

```
# Model A
model{
	for (j in 1:J1){
       		Y[j] ~ dbin(theta.adj[j], N[j])
		theta.adj[j] <- max(0.00001, min(theta, 0.99999))
	}
 
       	for (j in 1:J2){
       		Z[j] ~ dbern(p.adj[j])
		p.adj[j] <- max(0.00001, min(p[j], 0.99999))
		p[j] <- pbin(cut[j], theta, N[j+J1])  #Y<=cut
	}

	theta ~ dbeta(1,1)
}

# Model B

# Model C

# Model D

# Model E

# Model F

# Model G

```
