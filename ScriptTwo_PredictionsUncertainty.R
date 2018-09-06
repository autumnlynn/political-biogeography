##Load libraries
library(boot)
library(mgcv)
library(parallel)

##Required objects to be loaded in advance:
##1) eezs - list of EEZs
##2) mn.gam - fitted model object

##Calculate model predictions and uncertainty in model predictions following predict.gam help file
dnew <- data.frame(doy = 1:365, fid = sort(unique(data.mn$fid))[1], dum = 0) #Example covariate data for prediction from model without relative day of track term
#dnew <- data.frame(doy = 1:365, relday = 40:404, fid = sort(unique(data.mn$fid))[1], dum = 0) #Covariate data for prediction from model with relative day of track term; note that the correspondence between day of year and relative day of track is population/species-specific
p.hat <- predict(mn.gam, newdata = dnew, type = 'response', se.fit = TRUE) #Estimated probabilities
pred.link <- predict(mn.gam, newdata = dnew, type = 'link', se.fit = TRUE) #Estimated probabilities on link scale
Xp <- predict(mn.gam, dnew, type = 'lpmatrix') #Linear predictor matrix
sum.terms <- function(x){ #Function for summing linear predictor terms for each EEZ term group
  nc <- length(eezs) - 1 #Number of EEZ term groups (= number of EEZs - 1)
  nt <- length(x) / nc #Number of terms per EEZ term group
  res <- rep(NA, nc)
  for(i in 1:nc){
    res[i] <- sum(x[1:nt + (i - 1) * nt])
  }
  res
}
rmvn <- function(n, mu, sig){ #Function for generating multivariate normal random deviates of the model coefficients
  L <- mroot(sig)
  m <- ncol(L)
  t(mu + L %*% matrix(rnorm(m * n), m, n))
}
inv.link <- function(x){ #Function for calculating predicted probabilities (response) from predictions on link scale
  c(1, exp(x)) / (1 + sum(exp(x)))
}
postsamp <- function(x){ #Function for generating 'posterior sample' of predictions
  coef.r <- rmvn(1, coef(mn.gam), mn.gam$Vp) #Randomly sample coefficient vector
  prds <- t(apply(Xp, 1, '*', coef.r[1, ])) #Generate predictions on link scale for sample coefficient vector
  t(apply(prds, 1, sum.terms)) #Sum predictions by term group on link scale
}
nr <- 1e5 #Size of random sample for calculating confidence intervals
cl <- makeCluster(spec = 4, type = 'PSOCK') #Start parallel compute cluster
clusterExport(cl = cl, varlist = c('eezs', 'inv.link', 'mn.gam', 'nr', 'rmvn', 'sum.terms', 'Xp')) #Export required objects to parallel compute cluster
clusterCall(cl = cl, fun = library, package = 'mgcv', character.only = TRUE) #Load required libraries on parallel compute cluster
preds.link.r <- parLapply(cl = cl, X = 1:nr, fun = postsamp) #Generate 'posterior sample' of predictions on the link scale
na.rmv <- TRUE #Ignore bad individual link values and calculated probabilities?
preds.link.r <- array(unlist(preds.link.r), dim = c(dim(preds.link.r[[1]]), length(preds.link.r))) #Convert list to array
p.hat.r <- parApply(cl = cl, preds.link.r, c(1, 3), inv.link) #Convert predictions on link scale to probabilities
p.hat.r.lower50 <- apply(p.hat.r, c(1, 2), quantile, prob = 0.25, na.rm = na.rmv) #Lower bound of interquartile range of predicted probabilities across coefficient vector samples
p.hat.r.upper50 <- apply(p.hat.r, c(1, 2), quantile, prob = 0.75, na.rm = na.rmv) #Upper bound of interquartile range of predicted probabilities across coefficient vector samples
stopCluster(cl = cl) #Stop parallel compute cluster

##Calculate annual proportions of time in each EEZ
p.hat.ann <- apply(p.hat$fit, 2, function(x) sum(x) / dim(p.hat$fit)[1]) #Average probabilities across days for each EEZ, derived from model estimates
names(p.hat.ann) <- eezs
if(round(sum(p.hat.ann), 15) != 1) stop('p.hat.ann does not sum to 1')
p.hat.r.ann <- apply(p.hat.r, c(1, 3), function(x) sum(x) / dim(p.hat.r)[2]) #Average probabilities across days for each EEZ and replicate
idx <- which(apply(p.hat.r.ann, 2, function(x) all(!is.na(x)))) #Indexes of replicates with values for all EEZs
message(paste(round(length(idx) / dim(p.hat.r)[3] * 100, 0), '% of ', names(ann.probs)[wspc], ' replicates had values for all EEZs.', sep = ''))
p.hat.r.ann.lower50 <- apply(p.hat.r.ann[, idx], 1, quantile, prob = 0.25, na.rm = FALSE) #Lower bound of inter-quartile range of across-day average probabilities across replicates for each EEZ
p.hat.r.ann.upper50 <- apply(p.hat.r.ann[, idx], 1, quantile, prob = 0.75, na.rm = FALSE) #Upper bound of interquartile range of across-day average probabilities across replicates for each EEZ
