# Libraries ####
library(mgcv)
source("demons.R")

# Wrapper for standardizing tree ring data ####
# tra: the tree-ring array data structure containing the data to be analysed
# model: which effects (individual, time, age) to include in the model?
# form: are the effects added together or multiplied together
# error: is the data drawn from a normal additive PDF or a multiplicative log-normal PDF
# method: which algorithm should we use to standardize the data?
# sparse: use sparse list representation of data to reduce memory overhead

standardize_tra <- function(tra, model=list(I=FALSE, T=TRUE, A=TRUE), form="multiplicative", error="lnorm", method="sfs", sparse=TRUE, post_hoc=TRUE, ...)
{
  
  # Exception handling
  if (ifelse(is.data.frame(tra), sum(tra$G <= 0), sum(tra[tra<=0], na.rm=TRUE) > 0))
  {
    # Raise a warning if negative values found for multiplicative models
    if (form == "multiplicative")
    {
      warning("Zero or negative values found in data when using a multiplicative model. Estimated effects will not be stable.")
    }
    
    # Raise a warning if negative or 0 values found for lognormal errors
    if (error == "lnorm")
    {
      warning("Zero or negative values found in data when using a lognormal error form. Model fitting will fail, as non-positive values can not be generated using a log-normal distribution.")
    }
  }  
  
  # Raise a warning if error family and form do not match
  if (
      (form=="additive" & error=="lnorm")
      |
      (form=="multiplicative" & error=="norm")
  )
  {
    warning("Model form and error distribution are an unrealistic match. Be sure to check the model residuals. Model fitting problems may arise.")
  }

 
  if (method=="likelihood")
  {
    out <- standardize_likelihood(tra, model, form, error, sparse, ...)
  }
  else if(method == "least_squares")
  {
    out <- standardize_least_squares(tra, model, form, error, sparse, ...)
  }
  else if(method == "sfs")
  {
    #out <- standardize_sfs(tra, model, form, error, sparse, ...)
    out <- standardize_tsfs(tra, model, form, error, sparse, ...)
  }
  else if(method == "rcs")
  {
    out <- standardize_rcs(tra, model, form, error, sparse, ...)
  }
  else if(method == "gam")
  {
    out <- standardize_gam(tra, model, form, error, sparse, ...)
  }
  
  # Check for 3 effect model
  if (sum(unlist(model))==3)
  {
    if (post_hoc)
    {
      out$effects <- post_hoc_intercession(out$effects, out$tra, sparse, form)
      warning("Three effect model selected. Post-hoc selection was used to stabilize parameter estimates.")
    } else {
      warning("Three effect model selected. Parameter estimates are wildly unreliable. Consider using post-hoc selection.")
    }
  }
  
  
  return (out)
}

# Estimating and removing effects ####

# Naive estimate of a single effect (analagous to constructing regional curve or standardized chronology)
est_effect <- function (tra, factor.dim, mean_type="arithmetic", sparse=TRUE)
{  
  if(sparse)
  {
    id <- factor.dim + 1
    
    estimate_effect <- function(id_i)
    {
      data <- tra[tra[[id]]==id_i, "G"]
      if (mean_type=="geometric"){
        est_effect <- geomMean(data) 
      } else {
        est_effect <- mean(data, na.rm=TRUE)
      }
      return(est_effect)
    }
    
    est_effect <- sapply(levels(tra[[id]]), estimate_effect)
    
  }
  else
  {
    if (mean_type=="geometric"){
      est_effect <- apply(tra, factor.dim, geomMean) 
    } else {
      est_effect <- apply(tra, factor.dim, mean, na.rm=TRUE)
    }
  }
  
  return (est_effect)
}

# Remove an effect from a tree ring array
remove_effect <- function (tra, effect, factor.dim, form="multiplicative", sparse=TRUE)
{
  if(sparse)
  {
    id <- factor.dim + 1
        
    removed_tra <- tra
    
    for (effect_id in names(effect))
    {
      relevant_rows <- tra[[id]]==effect_id
      if (form=="additive")
      {
        removed_tra[relevant_rows,"G"] <- removed_tra[relevant_rows,"G"] - effect[effect_id]
      }
      else
      {
        removed_tra[relevant_rows,"G"] <- removed_tra[relevant_rows,"G"] / effect[effect_id]
      }
    }
  }
  else
  {
    FUN <- ifelse(form =="additive","-", "/")
    
    removed_tra <- sweep(tra, factor.dim, effect, FUN)
  }
  
  return (removed_tra)
  
}

# Add dummy effect vectors if some are missing
pad_effects <- function(effects, tra, form="multiplicative", sparse=TRUE)
{
  # Set the value to fill dummy coefficients with
  if (form=="multiplicative"){
    na.value <- 1
  } else {
    na.value <- 0
  }
  
  # Initialize dummy effects lists
  new_effects <- list(I=NA, T=NA, A=NA)
  
  if(sparse)
  {
    # Fill empty values  
    for(i in c("I", "T", "A")){
      
      j <- c(I="i", T="t", A="a")[i]
      
      if (length(effects[[i]] > 0)){
        new_effects[[i]] <-  effects[[i]]
        
        # Fill in dummy levels
        if(length(effects[[i]]) < nlevels(tra[[j]]))
        {
          missing_names <- levels(tra[[j]])[!(levels(tra[[j]]) %in% names(effects[[i]]))]
          missing_effects <- rep(na.value, times=length(missing_names))
          names (missing_effects) <- missing_names
          new_effects[[i]] <- c(effects[[i]], missing_effects)
        }
        
      } else {
        tra_dim <- which(c("I", "T", "A")==i)
        
        new_effects[[i]] <- rep.int(na.value, nlevels(tra[[tra_dim + 1]]))
        effect_names <- levels(tra[[tra_dim + 1]])
        if(i=="T" | i=="A")
        {
          effect_names <- as.numeric(as.character(effect_names))
        }
        names(new_effects[[i]]) <- sort(effect_names)
      }
    }
  }
  else
  {
    # Fill empty values  
    for(i in c("I", "T", "A")){
      
      j <- c(I=1, T=2, A=3)[i]
      
      if (length(effects[[i]] > 0)){
        new_effects[[i]] <-  effects[[i]]
        
        # Fill in dummy levels
        if(length(effects[[i]]) < dim(tra)[j])
        {
          missing_names <- dimnames(tra)[[j]][!(dimnames(tra)[[j]] %in% names(effects[[i]]))]
          missing_effects <- rep(na.value, times=length(missing_names))
          names (missing_effects) <- missing_names
          new_effects[[i]] <- c(effects[[i]], missing_effects)
        }
      } else {
        tra_dim <- which(c("I", "T", "A")==i)
        
        new_effects[[i]] <- rep.int(na.value, dim(tra)[tra_dim])
        names(new_effects[[i]]) <- dimnames(tra)[[tra_dim]]
      }
    }
  }
  
  return(new_effects)
}

# Correctly order effect vectors
sort_effects <- function(effects, tra, sparse=TRUE)
{
  if(sparse)
  {
    # Only ascending sorting as no "standard" order is retained
    sorted_effects <- list()
    
    # Sort I
    effect_names <- names(effects$I)
    sorted_effects$I <- effects$I[sort(effect_names)]
    
    # Sort T and A
    for (j in c("T", "A"))
    {
      effect_names <- as.numeric(names(effects[[j]]))
      sorted_effects[[j]] <- effects[[j]][as.character(sort(effect_names))]
    }
  }
  else
  {
    sorted_effects <- effects
    
    for (i in 1:3) 
    {
      correct_order <- dimnames(tra)[[i]]
      
      sorted_effects[[i]] <- effects[[i]][correct_order]
      
    }
  }
  
  return(sorted_effects)
}

# Rescale effect vectors to canonical form
rescale_effects <- function (effects, form="multiplicative")
{
  # Multiplicative models should be scaled such that the geometric mean of the secondary effects is 1
  # So log transform, set mean to 0, then unlog
  if (form=="multiplicative")
  {
    effects <- lapply(effects, log)
  }
  
  mean_effects <- lapply(effects, mean, na.rm=T)  
  
  # Scale I and T to mean of 0
  # Scale A so sum of effects stays the same
  
  # If A is missing, leave effects uncscaled
  if(!sum(!is.null(effects[[3]])))
  {
    return (effects)
  }
  
  # I
  if (sum (!is.null(effects[[1]]))){
    effects[[1]] <- effects[[1]]-mean_effects[[1]]
    effects[[3]] <- effects[[3]]+mean_effects[[1]]
  }
  
  # T
  if (sum (!is.null(effects[[2]]))){
    effects[[2]] <- effects[[2]]-mean_effects[[2]]
    effects[[3]] <- effects[[3]]+mean_effects[[2]]
  }
  
  if (form=="multiplicative")
  {
    effects <- lapply(effects, exp)
  }
  
  return (effects)
} 

# Model fit statistics ####

# Compute all the relevant model fit statistics for fixed-effects standardization
model_fit_tra <- function(effects, tra, model, form, error, sparse, method="sfs", k=NA)
{
  
  fit <- list()
  
  # Need to compute appropriate model fit stats for null models
  if (sum(unlist(model))==0)
  {
    # Predictions of the null model are the null value
    fit$predicted <- tra
    
    if (sparse) 
    {
      if (form=="additive")
      {
        fit$predicted$G <- 0
      } else
      {
        fit$predicted$G <- 1
      }
    } else 
    {
      if (form=="additive")
      {
        fit$predicted[!is.na(fit$predicted)] <- 0
      } else
      {
        fit$predicted[!is.na(fit$predicted)] <- 1
      }
    }
    
    # Residuals of the null model are the observed data
    fit$residuals <- tra
    
  } else {
    fit$predicted <- predicted_tra(effects, tra, form, sparse)
    fit$residuals <- residuals_tra (tra, fit$predicted, error, sparse)
  }
  
  fit$n <- n_tra(tra, sparse)
  if (method=="gam")
  {
    fit$k <- k
  } else{
    fit$k <- k_tra(tra, model, sparse)
  }
  fit$sigma <- sigma_tra(fit$residuals, error, sparse)
  
  fit$tss <- tss_tra(tra, error, sparse)
  fit$rss <- rss_tra(fit$residuals, error, sparse)
  fit$llh <- llh_tra(fit$residuals, error, sparse)
  
  fit$Rsq <- Rsq_tra(fit$rss, fit$tss)
  fit$adj.Rsq <- adj.Rsq_tra(fit$rss, fit$tss, fit$n, fit$k)
  
  fit$AIC <- AIC_tra(fit$llh, fit$k)
  fit$AICc <- AICc_tra(fit$llh, fit$k, fit$n)
  fit$BIC <- BIC_tra(fit$llh, fit$k, fit$n)
  
  return(fit)
}

# Predicted values
predicted_tra <- function (effects, tra, form, sparse)
{
  if(sparse)
  {
    predicted <- tra
    
    for (r in 1:nrow(predicted))
    {
      i <- as.character(tra[r, "i"])
      t <- as.character(tra[r, "t"])
      a <- as.character(tra[r, "a"])
      
      if (form=="additive")
      {
        predicted[r, "G"] <- effects$I[i] + effects$T[t] + effects$A[a]
      }
      else
      {
        predicted[r, "G"] <- effects$I[i] * effects$T[t] * effects$A[a]
      }
    }
    
  }
  else
  {
    predicted <- tra
    
    filled_cells <- which(!is.na(tra), arr.ind=TRUE)
    
    for (cell in 1:nrow(filled_cells))
    {
      i <- filled_cells[cell, 1]
      t <- filled_cells[cell, 2]
      a <- filled_cells[cell, 3]
      
      if (form=="additive")
      {
        predicted[i,t,a] <- effects$I[i] + effects$T[t] + effects$A[a]
      }
      else
      {
        predicted[i,t,a] <- effects$I[i] * effects$T[t] * effects$A[a]
      }
    }
  }
  
  
  return(predicted)
}

# Residuals
residuals_tra <- function (tra, predicted, error, sparse)
{
  if(sparse)
  {
    residuals <- tra
    
    if (error=="norm")
    {
      residuals$G <- tra$G - predicted$G
    }
    else
    {
      residuals$G <- tra$G / predicted$G
    }
  }
  else
  {
    if (error=="norm")
    {
      residuals <- tra - predicted
    }
    else
    {
      residuals <- tra / predicted
    }
  }
  
  return (residuals)
}

# Number of data points in the model
n_tra <- function (tra, sparse)
{
  if(sparse)
  {
    n <- nrow(tra)
  }
  else
  {
    n <- sum(!is.na(tra))
  }
  
  return(n)
}

# Number of parameters estimated (k)
k_tra <- function (tra, model, sparse)
{
  if(sparse)
  {
    # One parameter automatically for estimate of error
    k <- 1
    
    # One parameter is estimated for each index of the effect vectors
    if (model$I)
    {
      k <- k + nlevels(tra$i)
    }
    if (model$T)
    {
      k <- k + nlevels(tra$t)
    }
    if (model$A)
    {
      k <- k + nlevels(tra$a)
    }
    
    # Information about some parameters is lost due to rescaling (dummy variable trap)
    num_effects <- sum(unlist(model))
    
    k <- ifelse (num_effects > 0, k - (num_effects-1), 0)
  }
  else
  {
    # One parameter automatically for 
    k <- 1
    
    # One parameter is estimated for each index of the effect vectors
    if (model$I)
    {
      k <- k + dim(tra)[1]
    }
    if (model$T)
    {
      k <- k + dim(tra)[2]
    }
    if (model$A)
    {
      k <- k + dim(tra)[3]
    }
    
    # Information about some parameters is lost due to rescaling (dummy variable trap)
    num_effects <- sum(unlist(model))
    
    k <- ifelse (num_effects > 0, k - (num_effects-1), 0)
  }
  
  return(k)
}

# Calculate sigma, the level of dispersal in the noise PDF as the RMSE (the standard deviation of the residuals)
sigma_tra <- function (residuals, error, sparse)
{
  if(sparse)
  {
    val <- residuals$G
    if (error=="lnorm")
    {
      val <- log(val)
    }
    
    # sigma <- sqrt(mean(val^2))
    sigma <- sd(val)
  }
  else
  {
    val <- residuals[!is.na(residuals)]
    if (error=="lnorm")
    {
      val <- log(val)
    }
    
    # sigma <- sqrt(mean(val^2))
    sigma <- sd(val)
  }
  
  return(sigma)
}

# Total sum of squares
tss_tra <- function (tra, error, sparse)
{
  if(sparse)
  {
    val <- tra$G
    if (error=="lnorm")
    {
      val <- log(val)
    }
    
    mean_val <- mean(val)
    tss <- sum((val-mean_val)^2)
  }
  else
  {
    val <- tra[!is.na(tra)]
    if (error=="lnorm")
    {
      val <- log(val)
    }
    
    mean_val <- mean(val)
    tss <- sum((val-mean_val)^2)
  }
  
  return(tss)
}

# Residual sum of squares
rss_tra <- function (residuals, error, sparse)
{
  if(sparse)
  {
    val <- residuals$G
    if (error=="lnorm")
    {
      val <- log(val)
    }
    
    mean_val <- mean(val)
    rss <-sum((val-mean_val)^2)
  }
  else
  {
    val <- residuals[!is.na(residuals)]
    if (error=="lnorm")
    {
      val <- log(val)
    }
    
    mean_val <- mean(val)
    rss <- sum((val-mean_val)^2)
  }
  
  return(rss)
}

# Log-likelihood
llh_tra <- function (residuals, error, sparse)
{
  if(sparse)
  {
    sigma <- sigma_tra(residuals, error, sparse)
    
    val <- residuals$G
    
    # Likelihood is proportional to the probability of observing the data, given the parameters
    if(error=="norm")
    {
      llh <- sum(dnorm(val, sd=sigma, log=TRUE))
    }
    else
    {
      llh <- sum(dlnorm(val, sd=sigma, log=TRUE))
    }
  }
  else
  {
    sigma <- sigma_tra(residuals, error, sparse)
    
    val <- residuals[!is.na(residuals)]
    
    # Likelihood is proportional to the probability of observing the data, given the parameters
    if(error=="norm")
    {
      llh <- sum(dnorm(val, sd=sigma, log=TRUE))
    }
    else
    {
      llh <- sum(dlnorm(val, sd=sigma, log=TRUE))
    }
  }
  
  return(llh)
  
}

# R^2
Rsq_tra <- function (rss, tss)
{
  Rsq <- 1-rss/tss
  return(Rsq)
}

# Adjusted R^2
adj.Rsq_tra <- function (rss, tss, n, k)
{
  # p doesn't include estimates of the variance
  p <- k - 1
  adj.Rsq <- 1-(rss/tss)*(n-1)/(n-p-1)
  return(adj.Rsq)
}

# AIC
AIC_tra <- function (llh, k)
{
  AIC <- -2*llh + 2*k
  return(AIC)
}

# AICc
AICc_tra <- function (llh, k, n)
{
  AICc <- -2*llh + 2*k + 2*k*(k+1)/(n-k-1)
  return(AICc)
}

# BIC
BIC_tra <- function (llh, k, n)
{
  BIC <- -2*llh + k*log(n)
  return(BIC)
}

# Regional curve standardization ####
# effect_order: the order in which effects are sequentially estimated
standardize_rcs <- function(tra, model=list(I=FALSE, A=TRUE, T=TRUE), form="multiplicative", error="lnorm", sparse=TRUE)
{
  
  # Convert the tree-ring array to the appropriate form (sparse/full)
  if (sparse)
  {
    if (!is.data.frame(tra))
    {
      tra <- sparse_tra(tra)
    }
  } else
  {
    if (is.data.frame(tra))
    {
      tra <- unsparse_tra(tra)
    }
  }
  
  # Select appropriate type of mean
  if (error=="lnorm"){
    mean_type <- "geometric"
  } else {
    mean_type <- "arithmetic"
  }
  
  # Determine effect order from order in which I, T, A is listed
  inc_effects <- names(model[unlist(model)==TRUE])
  name_dim_dict <- c(I="I", T="T", A="A")
  effect_order <- match(inc_effects, name_dim_dict)
  
  # Make a dummy list of effects
  effects <- list(I=NULL, T=NULL, A=NULL)
  
  # Make a dummy tree-ring array
  working_tra <- tra
  
  # Estimate the effects one at a time
  for (effect in effect_order){
    # Estimate an effect    
    effects[[effect]] <- est_effect(working_tra, effect, mean_type, sparse)
    
    # Remove the effect
    working_tra <- remove_effect(working_tra, effects[[effect]], effect, form, sparse)
  }
  
  # Fill in dummy values for effects not estimated
  effects <- pad_effects(effects, tra, form, sparse)
  
  # Make sure effects are in the right order
  effects <- sort_effects(effects, tra, sparse)
  
  # Rescale the effects to standard form
  effects <- rescale_effects(effects, form)
  
  # Compute model fit statistics
  fit <- model_fit_tra (effects, tra, model, form, error, sparse)
    
  # Record model fitting settings
  settings <- list(model=model, form=form, error=error, sparse=sparse, method="rcs")
  
  out <- list(effects=effects, tra=tra, fit=fit, settings=settings)
  
  return (out)   
}

# Signal-free regional curve standardization ####
standardize_sfs <-function (tra, model=list(I=FALSE, T=TRUE, A=TRUE), form="multiplicative", error="lnorm", sparse=TRUE, cor_threshold=0.999999)
{
  # Convert the tree-ring array to the appropriate form (sparse/full)
  if (sparse)
  {
    if (!is.data.frame(tra))
    {
      tra <- sparse_tra(tra)
    }
  } else
  {
    if (is.data.frame(tra))
    {
      tra <- unsparse_tra(tra)
    }
  }
  
  # Select appropriate type of mean
  if (error=="lnorm"){
    mean_type <- "geometric"
  } else {
    mean_type <- "arithmetic"
  }
  
  # Determine effect order from order in which I, T, A is listed
  # Effect order needs to be of length 2 to work  
  inc_effects <- names(model[unlist(model)==TRUE])
  name_dim_dict <- c(I="I", T="T", A="A")
  effect_order <- match(inc_effects, name_dim_dict)
  
  # Loop controls
  converged <- FALSE
  iteration <- 0
  
  working_tra <- tra
  
  while (!converged){
    
    # Upkeep counters
    last_tra <- working_tra
    iteration <- iteration +1
    print (paste("Iteration", iteration))
    
    # Estimate the values for the first dimension
    est_1 <- est_effect(working_tra, effect_order[1], mean_type, sparse)
    #print(est_1)
    
    # Remove those effects temporarily
    intermediate_tra <- remove_effect (working_tra, est_1, effect_order[1], form, sparse)
    
    # Estimate values for the second dimensions
    est_2 <- est_effect(intermediate_tra, effect_order[2], mean_type, sparse)
    #print(est_2)
    
    # Remove those effects from the working data
    working_tra <- remove_effect (working_tra, est_2, effect_order[2], form, sparse)
    
    # Check for convergence. Use the log-correlation if the error term is suspected to be multiplicative lognormal
    if (error=="norm"){
      if (sparse)
      {
        conv_cor <- cor(working_tra$G, last_tra$G) 
      }else{
        conv_cor <- cor(working_tra, last_tra,  "complete.obs") 
      }
      
    } else {
      if (sparse)
      {
        conv_cor <- cor(log(working_tra$G), log(last_tra$G))
      }else{
        conv_cor <- cor(log(working_tra), log(last_tra),  "complete.obs")
      }
    }
        
    if (conv_cor>=cor_threshold){
      converged <- TRUE
    }
    
  }
  
  # Create storage for the estimated effects
  effects <- vector(mode="list", length=3)
  
  # Dummy missing effects
  all.dim <- c(1,2,3)
  
  miss.dim <- all.dim[-intersect(effect_order, all.dim)]
  miss_effect <- rep.int(1, dim(tra)[[miss.dim]])
  names (miss_effect) <- dimnames(tra)[[miss.dim]]
  
  # Primary chronology is mean of converged working TRA
  prim_effect <- est_effect(working_tra, effect_order[1], mean_type, sparse)
  
  # Secondary chronology is mean of original data with primary effects removed
  sec.series <- remove_effect (tra, prim_effect, effect_order[1], form, sparse)
  sec_effect<- est_effect(sec.series, effect_order[2], mean_type, sparse)
  
  # Compile the effects in the approriate order
  effects <- list(prim_effect, sec_effect, miss_effect)
  
  effect_order <- order(c(effect_order[1], effect_order[2], miss.dim))
  effects <- effects[effect_order]
  
  # Make sure effects are in the right order
  effects <- sort_effects(effects, tra, sparse)
  
  # Rescale the effects  to standard form
  effects <- rescale_effects(effects, form)
  
  # Compute model fit statistics
  fit <- model_fit_tra (effects, tra, model, form, error, sparse)
  
  # Record model fitting settings
  settings <- list(model=model, form=form, error=error, sparse=sparse, method="old_sfs")
  
  out <- list(effects=effects, tra=tra, fit=fit, settings=settings)
  
  return (out)
}

# Truly signal-free regional curve standardization ####
# Cleans up SF-RCS algorithm and allows expansion to N dimensions

standardize_tsfs <- function (tra, model=list(I=FALSE, T=TRUE, A=TRUE), form="multiplicative", error="lnorm", sparse=TRUE, cor_threshold=0.999999)
{
  
  # Convert the tree-ring array to the appropriate form (sparse/full)
  if (sparse) 
  {
    if (!is.data.frame(tra))
    {
      tra <- sparse_tra(tra)
    }
  } else
  {
    if (is.data.frame(tra))
    {
      tra <- unsparse_tra(tra)
    }
  }
  
  # Select appropriate type of mean
  if (error=="lnorm"){
    mean_type <- "geometric"
  } else {
    mean_type <- "arithmetic"
  }
  
  # Create storage for the estimated effects
  effects <-vector(mode="list", length=3)
  names (effects) <- c("I","T","A")
  
  # Dummy starting effects
  for (i in 1:3){
    if(sparse)
    {
      dim_i <- nlevels(tra[[i+1]])
    } else
    {
      dim_i <- dim(tra)[[i]]
    }
    
    if (form=="additive")
    {
      effects[[i]] <-  rep.int(0,  dim_i)
    } else
    {
      effects[[i]] <-  rep.int(1,  dim_i)
    }
    
    if(sparse)
    {
      names(effects[[i]]) <- levels(tra[[i+1]])
    } else
    {
      names(effects[[i]]) <- dimnames(tra)[[i]]
    }
  }
  
  # Determine effect order from order in which I, T, A is listed
  inc_effects <- names(model[unlist(model)==TRUE])
  name_dim_dict <- c(I="I", T="T", A="A")
  effect_order <- match(inc_effects, name_dim_dict)
  
  # Loop controls
  
  # Don't attempt to fit null models
  if (length(inc_effects)==0)
  {
    converged <- TRUE 
  } else 
  {
    converged <- FALSE
  }
  iteration <- 0
  working_tra <- tra
  
  while (!converged){
    
    # Upkeep counters
    last_tra <- working_tra
    iteration <- iteration +1
    print (paste("Iteration", iteration))
    
    for (j in effect_order){
      
      # Estimate the effects across each dimension
      est_j <- est_effect(working_tra, j, mean_type, sparse)
      
      # Remove them from the signal-free data
      working_tra <- remove_effect (working_tra, est_j, j, form, sparse)
      
      # Combine them with previously determined effects for that dimension
      if (form == "additive")
      {
        effects[[j]] <- effects[[j]] + est_j
      } else
      {
        effects[[j]] <- effects[[j]] * est_j
      }
    }
    
    # Check for convergence. Use the log-correlation if the error term is suspected to be multiplicative lognormal    
    if (error=="norm"){
      if (sparse)
      {
        conv_cor <- cor(working_tra$G, last_tra$G) 
      }else{
        conv_cor <- cor(working_tra, last_tra,  "complete.obs") 
      }
      
    } else {
      if (sparse)
      {
        conv_cor <- cor(log(working_tra$G), log(last_tra$G))
      }else{
        conv_cor <- cor(log(working_tra), log(last_tra),  "complete.obs")
      }
    }
    
    if (conv_cor>=cor_threshold){
      converged <- TRUE
    }
  
  }
  
  # Make sure effects are in the right order
  effects <- sort_effects(effects, tra, sparse)
  
  # Rescale the effects to standard form
  effects <- rescale_effects(effects, form)

  # Compute model fit statistics
  fit <- model_fit_tra (effects, tra, model, form, error, sparse)
  
  # Record model fitting settings
  settings <- list(model=model, form=form, error=error, sparse=sparse, method="sfs")
  
  out <- list(effects=effects, tra=tra, fit=fit, settings=settings)
  
  return (out)
}

# Maximum likelihood fixed effects standardization ####
# Fits models using maximum likelihood
# Searches for solutions with simulated annealing

# Flattening effects to a vector for use with optim()
flatten_effects <- function(effects, model)
{
  flat_effects <- vector()
  
  if (model$I)
  {
    effects_I <- dummy_effects$I
    names(effects_I) <- paste("I", names(dummy_effects$I), sep="@")
    flat_effects <- c(flat_effects, effects_I)
  }
  if (model$T)
  {
    effects_T <- dummy_effects$T
    names(effects_T) <- paste("T", names(dummy_effects$T), sep="@")
    flat_effects <- c(flat_effects, effects_T)
  }
  if (model$A)
  {
    effects_A <- dummy_effects$A
    names(effects_A) <- paste("A", names(dummy_effects$A), sep="@")
    flat_effects <- c(flat_effects, effects_A)
  }
  return(flat_effects)
}

unflatten_effects <- function (flat_effects)
{
  effects <- list()
  
  id_code <- substr(names(flat_effects), 1, 1)
  names(flat_effects) <- substr(names(flat_effects), 3, length(names(flat_effects)))
  
  effects$I <- flat_effects[id_code=="I"]
  effects$T <- flat_effects[id_code=="T"]
  effects$A <- flat_effects[id_code=="A"]
  
  return(effects)
}

standardize_mle <- function(tra, model=list(I=FALSE, T=TRUE, A=TRUE), form="multiplicative", error="lnorm", sparse=TRUE, ...)
{
  # Convert the tree-ring array to the appropriate form (sparse/full)
  if (sparse) 
  {
    if (!is.data.frame(tra))
    {
      tra <- sparse_tra(tra)
    }
  } else
  {
    if (is.data.frame(tra))
    {
      tra <- unsparse_tra(tra)
    }
  }
  
  # Determine starting estimates for the effects
  dummy_effects <- pad_effects(list(), tra, form, sparse)
  flat_effects0 <- flatten_effects(dummy_effects, model)
  
  # Likelihood function, optimize this!
  fes_likelihood <- function(flat_effects)
  {
    # Unflatten effects
    effects <- unflatten_effects(flat_effects)
    
    # Pad the effects
    effects <- pad_effects(effects, tra, form, sparse)
    
    # Find the predicted values
    predicted <- predicted_tra(effects, tra, form, sparse)
      
    # Find the residuals
    residuals <- residuals_tra(tra, predicted, error, sparse)
    
    # Find the likelihood
    llh <- llh_tra(residuals, error, sparse)
      
    # Return the negative likelihood (optim is a minimizer)
    return(-llh)
  }
  
  # Optimize the model
  mle_solution <- optim(flat_effects0, fes_likelihood, ...)
  
  # Report the optimizer's output
  print(mle_solution[2:5])
  
  # Extract the effects
  effects <- unflatten_effects(mle_solution$par)
  
  # Fill in dummy values for effects not estimated
  effects <- pad_effects(effects, tra, form, sparse)
  
  # Make sure effects are in the right order
  effects <- sort_effects(effects, tra, sparse)
    
  # Rescale the effects to standard form
  effects <- rescale_effects(effects, form)
  
  # Compute model fit statistics
  fit <- model_fit_tra (effects, tra, model, form, error, sparse)
  
  # Record model fitting settings
  settings <- list(model=model, form=form, error=error, sparse=sparse, method="likelihood", ...)
  
  out <- list(effects=effects, tra=tra, fit=fit, settings=settings)
  
  return(out)
}

# GAM fixed effects standardization ####

# Main gam function
standardize_gam <- function (tra, model=list(I=FALSE, T=TRUE, A=TRUE), form="multiplicative", error="lnorm", sparse=TRUE, age_k=10, ...)
{
  # Confirm that form and error match, otherwise GAMs can't be used
  if (
    (form=="additive" & error=="lnorm")
    |
      (form=="multiplicative" & error=="norm")
  )
  {
    simpleError("Model form and error distribution are an unrealistic match. Generalized additive models cannot be used for this configuration.")
  }
  
  # Convert the tree-ring array to the appropriate form (sparse/full)
  if (sparse) 
  {
    if (!is.data.frame(tra))
    {
      tra <- sparse_tra(tra)
    }
  } else
  {
    # Regression uses table form data
    simpleError("GAM standardization can only handle sparse tree-ring arrays.")
  }
    
  # Construct formula for regression
  growth_formula <- as.formula(make_gam_formula(model))
  
  # max_k determines the maximum flexibility of the spline fitting the age effect
  # Increased max_k greatly increases computation time
  # Absolute maximum flexibility is:
  # max_k <- nlevels(tra$a)
  
  print ("Model initialized.")
  print (growth_formula)
  
  # If form is additive and error is additive and normal:
  # Use GAM regression with a log-link
  if (form == "additive")
  {
    family <- gaussian(link="identity")
  }
  
  # If form is multiplicative and error is multiplicative log-normal:  
  # Use GAM regression with a log-link
  if (form == "multiplicative")
  {
    # Extremely slow, but more correct
    #family <- Gamma(link="log")
    # Log transformed response variable hack
    # Breaks gam-estimated AIC
    tra$G <- log(tra$G)
    family <- gaussian(link="identity")
  }
  
  # Estimate the growth model
  growth_model <- gam (growth_formula, family=family, data=tra, ...)
  print ("Generalized additive model used to standardize data.")
  
  # Extract estimates of the effect
  # Correct way
  #effects <- extract_effects_gam(growth_model, model, form, tra)
  
  # Log-transformation hack
  effects <- extract_effects_gam(growth_model, model, form="additive", tra)
  if (form=="multiplicative")
  {
    effects <- lapply(effects, exp)
    tra$G <- exp(tra$G)
  }
  print ("Effects extracted.")
  
  # Fill in dummy values for effects not estimated
  effects <- pad_effects(effects, tra, form, sparse)
  
  # Make sure effects are in the right order
  effects <- sort_effects(effects, tra, sparse)
  
  # Rescale the effects to standard form
  effects <- rescale_effects(effects, form)
  
  # Compute model fit statistics
  # k, the number of parameters fit, must be gathered from the model fitting function
  # k is ~= the effective degrees of freedom + 1
  k <- sum(growth_model$edf) + 1
  fit <- model_fit_tra (effects, tra, model, form, error, sparse, method="gam", k=k)
  print("Fit computed.")
  
  # Record model fitting settings
  settings <- list(model=model, form=form, error=error, sparse=sparse, method="gam", ...) #max_k=max_k, ...)
  
  out <- list(effects=effects, tra=tra, fit=fit, settings=settings)
  
  return (out)
  
}

# Formula construction for GAM standardization
make_gam_formula <- function (model){
  dep.str <- "G"
  ind.str <- "0"
  
  if(model$I){
    ind.str <- paste(ind.str, "i", sep="+")
  }
  if(model$T){
    ind.str <- paste(ind.str, "t", sep="+")
  }
  if(model$A){
    # Change smoothing terms here
    ind.str <- paste(ind.str, "s(as.numeric(as.character(a)), k=age_k, ...)", sep="+")
  }
  
  # Combine the two sides of the formula
  formula.str <- paste(dep.str, ind.str, sep="~")
  
  return (formula.str)
}

# Extracting effects for gam models
extract_effects_gam <- function(growth_model, model, form, tra)
{
  
  # Extract coefficients for I and T directly
  raw_effects <- growth_model$coefficients
  
  effects <- list()
  effects$I <- raw_effects[substr(names(raw_effects), 1, 1)=="i"]
  effects$T <- raw_effects[substr(names(raw_effects), 1, 1)=="t"]

  names(effects$I) <- substr(names(effects$I), 2, length(names(effects$I)))
  names(effects$T) <- substr(names(effects$T), 2, length(names(effects$T)))
  
  # Find A by process of elimination
  # Generate predicted values
  dummy_data <- data.frame(i=levels(tra$i)[2], t=levels(tra$t)[2], a=levels(tra$a))
  predicted_by_age <- predict(growth_model, dummy_data)
  
  # Remove the known effects of time and individuals
  # Convoluted partial matching code to handle variable name truncation
  if (form=="additive")
  {
    base_line <- 0
    if (model$I)
    {
      base_line <- base_line + effects$I[which(!is.na(pmatch(names(effects$I), levels(tra$i)[2])))]
    }
    if (model$T)
    {
      base_line <- base_line + effects$T[which(!is.na(pmatch(names(effects$T), levels(tra$t)[2])))]
    }
    effects$A <- predicted_by_age - base_line
  } else 
  {
    base_line <- 1
    if (model$I)
    {
      base_line <- base_line * effects$I[which(!is.na(pmatch(names(effects$I), levels(tra$i)[2])))]
    }
    if (model$T)
    {
      base_line <- base_line * effects$T[which(!is.na(pmatch(names(effects$T), levels(tra$t)[2])))]
    }
    effects$A <- predicted_by_age / base_line
  }
  
  names(effects$A) <- levels(tra$a)
  
  return(effects)
}

