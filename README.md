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

## Illustrative Examples:

### 1. Survival data

``Data``

``JAGS Model``

### 2. Binomial data

``Data``

``JAGS Model``
