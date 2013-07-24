# Load dependencies ####
source ("coreTRA.R")
source("FRS.R")
source("RCS.R")
source("syntheticData.R")
source ("plotting.R")
source("multicollinearity.R")

library(pastecs)
library(fBasics)
library (bbmle)
library(reshape2)

# Jack pine site functions ####

getSite <- function (I.names){
  sites <- sapply(I.names, substr, 1, 3)
  
  df <- data.frame (Tree=I.names, Site=sites)
  
  return (df)
}

firstYear <- function (index, tra){
  full <- which(!is.na(tra[index,,]), arr.ind=T)
  first <- min (full[[1]])
  return (first)
}

getSiteStart <- function (tra){
  Q.names <- dimnames(tra)[[1]]
  
  # Get site code for each core
  sites <- getSite(Q.names)
  
  # Find first year of each core
  firstYears <- vector ()
  for (i in Q.names){
    firstYears[i] <- firstYear (i, tra)
  }
  
  siteStart <- vector()
  
  # Find first year for each site
  for (s in levels(sites$Site)){
    allTreesFirst <- firstYears [sites [which(sites$Site==s), "Tree"]]
    
    siteStart[s] <- min (allTreesFirst)
    
    # Convert from index to real year
    siteStart[s] <- dimnames(tra)[[2]][as.numeric(siteStart[s])]
  }
  
  return (siteStart)
}

evenAgedTRA <- function (rwl){
  simpleTRA <- rwl.to.tra(rwl)
  
  siteStarts <- getSiteStart (simpleTRA)
  
  siteDF <- getSite (dimnames(simpleTRA)[[1]])
  
  newTreeStarts <- vector()
  
  for (i in dimnames(simpleTRA)[[1]]){
    site <- siteDF[which (siteDF$Tree==i), "Site"]
    newTreeStarts[i] <- siteStarts[site]
  }
  
  evenTRA <- rwl.to.tra (rwl, newTreeStarts)
  
  return (evenTRA)
}

# Jack pine analysis ##################################################

jp <- read.rwl ("C:/Users/Homeuser/Desktop/Documents/Research/TreeRingArrays/Data/masterclean.rwl")

# Remove the year 2012 from the rwl table. Should not exist
jp <- jp[-which(rownames(jp)=="2012"),]

# Even aged
# Convert to TRA format
jp.tra <- evenAgedTRA (jp)
# jp.tra <- rwl.to.tra(jp)

# Convert 0 values to NA
jp.tra[which(jp.tra==0, arr.ind=T)] <- NA

# Report sample depth
print (plotSampleDepth(jp.tra, 1))
print (plotSampleDepth(jp.tra, 2))
print (plotSampleDepth(jp.tra, 3))

# Fit the data using the two alternate models and regression functions ####
# Optimal spans: QFA=0.0606, FA=0.0614
jp.gam.QFA <- standardize_frs(jp.tra, T,T,T, model_type="gam", span="auto")
jp.gam.FA <-  standardize_frs(jp.tra, F,T,T, model_type="gam", span="auto")

# Compile CV ####
Q_jp <- data.frame("QFA"=jp.gam.QFA$cv[[1]])
F_jp <- data.frame("QFA"=jp.gam.QFA$cv[[2]], "FA"=jp.gam.FA$cv[[2]])
A_jp <- as.data.frame(t(smartbind("QFA"=jp.gam.QFA$cv[[3]], "FA"=jp.gam.FA$cv[[3]])))


# Generate an age list
ageList <- rep.int(NA, nrow(Q_jp))
names(ageList) <- rownames (Q_jp)

for (i in 1:nrow(Q_jp)){
  ageList[i] <- max(which(!is.na(jp.tra[i,,]), arr.ind=T)[,2])
}

# Plotting ####

# General pdf saving parameters
# pdf_width=11
# pdf_height=8.5

# Colours: 
# GAM: Dark
# LM: Light
# 3-factor: Red
# 2-factor: Blue
# Uncorrected: Green

# Root save path
main_wd <- getwd()
figure_path <- paste(getwd(), "Figures", "Pinus_banksiana", sep="/")

# Saving path
# save_path <- paste(figure_path, "Ring width", sep="/")
# setwd(save_path)

# Plot Q
# Plot Q vs age
# pdf("Q_age.pdf", width=pdf_width, height=pdf_height)
Q_vs_age_plot (Q_jp,jp.tra) + scale_color_manual(values=c("firebrick4"))
# dev.off()

# Plot mean Q by year
# pdf("Q_year.pdf", width=pdf_width, height=pdf_height)
meanQ_by_year_plot (Q_jp,jp.tra)  + scale_color_manual(values=c("firebrick4"))
# dev.off()

# Plot F
# pdf("F.pdf", width=pdf_width, height=pdf_height)
plot.cv(F_jp, factor.type="F") + scale_color_manual(values=c("blue4", "firebrick4"))
# dev.off()

# Plot A
# pdf("A.pdf", width=pdf_width, height=pdf_height)
plot.cv(A_jp, factor.type="A")  + scale_color_manual(values=c("blue4", "firebrick4"))+ylab("Ring width (mm)")
# dev.off()

# Reset working directory
# setwd(main_wd)

# Diagnostics ####

# # Saving path
# save_path <- paste(figure_path, "Residuals", sep="/")
# setwd(save_path)

# Check if Q is correlated with age
my_cor.test<- function(...){
  return(cor.test(..., method="spearman"))
}

sapply(Q_jp, cor.test, y=ageList, method="spearman")

# Test if there is a trend in forcing or age
sapply(F_jp, my_cor.test, y=as.numeric(rownames(F_jp)))
sapply(A_jp, my_cor.test, y=as.numeric(rownames(A_jp)))

# # Check the residuals
# raw_jp <- data.frame(x=jp.tra[!is.na(jp.tra)])
# raw_jp [raw_jp==0] <- NA
# log_jp <-  log(raw_jp)
# resid_gam_jp_2 <- data.frame(x=residuals (jp.gam.FA$model))
# resid_gam_jp_3 <- data.frame(x=residuals (jp.gam.QFA$model))
# resid_lm_jp_2 <- data.frame(x=residuals (jp.lm.FA$model))
# resid_lm_jp_3 <- data.frame(x=residuals (jp.lm.QFA$model))
# 
# # Residual standard deviation
# sapply(log_jp, sd)
# sapply(resid_gam_jp_2, sd)
# sapply(resid_gam_jp_3, sd)
# sapply(resid_lm_jp_2, sd)
# sapply(resid_lm_jp_3, sd)
# 
# # Q-Q plots
# qq_raw <- as.data.frame(qqnorm(raw_jp$x))
# qq_log <- as.data.frame(qqnorm(log_jp$x))
# qq_gam_resid_2 <- as.data.frame(qqnorm(resid_gam_jp_2$x))
# qq_gam_resid_3 <- as.data.frame(qqnorm(resid_gam_jp_3$x))
# qq_lm_resid_2 <- as.data.frame(qqnorm(resid_lm_jp_2$x))
# qq_lm_resid_3 <- as.data.frame(qqnorm(resid_lm_jp_3$x))
# 
# pdf("qq_raw.pdf", width=pdf_width, height=pdf_height)
# ggplot(data=qq_raw, aes(x=x, y=y))+geom_point()+theme_bw()+ylab("Sample Quantiles")+xlab("Theoretical Quantiles")+ geom_abline(intercept = 
#                                                                                                                                  mean(qq_raw$y), slope = sd(qq_raw$y)) 
# dev.off()
# 
# pdf("qq_log.pdf", width=pdf_width, height=pdf_height)
# ggplot(data=qq_log, aes(x=x, y=y))+geom_point()+theme_bw()+ylab("Sample Quantiles")+xlab("Theoretical Quantiles")+ geom_abline(intercept = mean(qq_log$y), slope = sd(qq_log$y))
# dev.off()
# 
# pdf("qq_gam_2.pdf", width=pdf_width, height=pdf_height)
# ggplot(data=qq_gam_resid_2, aes(x=x, y=y))+geom_point()+theme_bw()+ylab("Sample Quantiles")+xlab("Theoretical Quantiles")+ geom_abline(intercept = mean(qq_gam_resid_2$y), slope = sd(qq_gam_resid_2$y))
# dev.off()
# 
# pdf("qq_gam_3.pdf", width=pdf_width, height=pdf_height)
# ggplot(data=qq_gam_resid_3, aes(x=x, y=y))+geom_point()+theme_bw()+ylab("Sample Quantiles")+xlab("Theoretical Quantiles")+ geom_abline(intercept = mean(qq_gam_resid_3$y), slope = sd(qq_gam_resid_3$y))
# dev.off()
# 
# pdf("qq_lm_2.pdf", width=pdf_width, height=pdf_height)
# ggplot(data=qq_lm_resid_2, aes(x=x, y=y))+geom_point()+theme_bw()+ylab("Sample Quantiles")+xlab("Theoretical Quantiles")+ geom_abline(intercept = mean(qq_lm_resid_2$y), slope = sd(qq_lm_resid_2$y))
# dev.off()
# 
# pdf("qq_lm_3.pdf", width=pdf_width, height=pdf_height)
# ggplot(data=qq_lm_resid_3, aes(x=x, y=y))+geom_point()+theme_bw()+ylab("Sample Quantiles")+xlab("Theoretical Quantiles")+ geom_abline(intercept = mean(qq_lm_resid_3$y), slope = sd(qq_lm_resid_3$y))
# dev.off()
# 
# # Density
# pdf("density_raw.pdf", width=pdf_width, height=pdf_height)
# ggplot(data=raw_jp, aes(x=x)) + geom_density(alpha=0.3, fill="black")+theme_bw()+xlab("Basal area increment (mm2)")+ylab("Density")
# dev.off()
# 
# foo <- data.frame ("Raw"=log_jp, "GAM-2"=resid_gam_jp_2, "GAM-3"=resid_gam_jp_3, "LM-2"=resid_lm_jp_2, "LM-3"=resid_lm_jp_3)
# names (foo) <- c("Raw", "GAM-2", "GAM-3", "LM-2", "LM-3")
# foo <- data.frame(x=unlist(foo), Model=rep(names(foo), each=nrow(foo)))
# 
# pdf("density_all.pdf", width=pdf_width, height=pdf_height)
# ggplot(data=foo, aes(x=x, fill=Model)) + geom_density(alpha=0.3)+theme_bw()+xlab("Log Basal area increment (mm2)")+geom_vline(x=0)+ylab("Density") + scale_fill_manual(values=c("blue4", "firebrick4", "deepskyblue",  "indianred1", "black"))
# dev.off()
# 
# bar <- data.frame ("Raw"=log_jp, "GAM-2"=resid_gam_jp_2, "GAM-3"=resid_gam_jp_3)
# names (bar) <- c("Raw", "GAM-2", "GAM-3")
# bar <- data.frame(x=unlist(bar), Model=rep(names(bar), each=nrow(bar)))
# 
# pdf("density_gam.pdf", width=pdf_width, height=pdf_height)
# ggplot(data=bar, aes(x=x, fill=Model)) + geom_density(alpha=0.3)+theme_bw()+xlab("Log Basal area increment (mm2)")+geom_vline(x=0)+ylab("Density") + scale_fill_manual(values=c("firebrick4", "deepskyblue", "black"))
# dev.off()

# Reset working directory
# setwd(main_wd)

# Carbon-13 analysis ####
raw_deltaC_Wood <- read.csv("C:/Users/Homeuser/Desktop/Documents/Research/TreeRingArrays/Data/tree_13c.csv")
raw_deltaC_Wood <- c13[1:3]

Ca <- read.csv("C:/Users/Homeuser/Desktop/Documents/Research/TreeRingArrays/Data/co2_record.csv")
Ca <- Ca[-3]
rownames(Ca) <- Ca[[1]]
Ca <- Ca[-1]
Ca <- data.frame(Ca=Ca[as.numeric(rownames(Ca))<=2010 & as.numeric(rownames(Ca))>=1877,])
rownames(Ca) <- 1877:2010

deltaC_Atm <- read.csv("C:/Users/Homeuser/Desktop/Documents/Research/TreeRingArrays/Data/13C_record.csv")
rownames(deltaC_Atm) <- deltaC_Atm[[1]]
deltaC_Atm <- deltaC_Atm[-1]
deltaC_Atm <- data.frame(DC=deltaC_Atm[as.numeric(rownames(deltaC_Atm))<=2010 & as.numeric(rownames(deltaC_Atm))>=1877,])
rownames(deltaC_Atm) <- 1877:2010

# # Convert format to a tree ring table
deltaC_Wood <-  dcast(raw_deltaC_Wood, Year~Site, mean)
deltaC_Wood[sapply(deltaC_Wood, is.nan)] <- NA
rownames(deltaC_Wood) <- deltaC_Wood[[1]]
deltaC_Wood <- deltaC_Wood[-1]

# Convert to DeltaC
DeltaC <- (deltaC_Atm[[1]]-deltaC_Wood)/(1000+deltaC_Wood)*1000

# Convert to iWUE
a <- 4.4
b <- 27
Wi <- Ca[[1]]/1.6*(1-(DeltaC-a)/(b-a))

# Run Factor Regression Standardization
# Span of 0.135
Wi_frs <- standardize_frs(rwl.to.tra(Wi), model_type="gam")

# Rescale CV so F is dominant
Q_Wi <- Wi_frs$cv[[1]]
F_Wi <- Wi_frs$cv[[2]]*geomMean(Wi_frs$cv[[3]])
A_Wi <- Wi_frs$cv[[3]]/geomMean(Wi_frs$cv[[3]])

# Find climate-growth links ####
# Preprocess forcing info from tree rings
forcing_QFA <- F_jp[1][as.numeric(rownames(F_jp))<=2007 & as.numeric(rownames(F_jp))>=1901,]
forcing_FA <- F_jp[2][as.numeric(rownames(F_jp))<=2007 & as.numeric(rownames(F_jp))>=1901,]
names(forcing_QFA) <- 1901:2007
names(forcing_FA) <- 1901:2007

# Load climate data
jp.clim <- read.csv("C:/Users/Homeuser/Desktop/Documents/Research/TreeRingArrays/Data/IgnaceClimate.csv")
foo <- melt(jp.clim, id.vars=c("Month", "Year"))
jp.clim.wide <- recast(foo, Year~variable+Month) 
full_clim <- as.data.frame(read.csv("C:/Users/Homeuser/Desktop/Documents/Research/TreeRingArrays/Data/full_clim.csv", header=T))
rownames(full_clim) <- full_clim[[1]]
full_clim <- full_clim[-1:-5]

# Preliminary investigations into good candidate  variables
detach (package:gam)
library(mgcv)

clim_fit <- vector()
for (var in 1:ncol(full_clim)){
  clim_var <- full_clim[[var]]
  try(prelim_gam <- gam(forcing_QFA~s(clim_var)))
  clim_fit <- c(clim_fit, summary(prelim_gam)$r.sq)
}
names(clim_fit) <- names(full_clim)
clim_fit <- sort(clim_fit, decreasing=T)

# Prelim variables, R^2 > 0.1:
# 1, 9, 10, 11, 12, 15, 19, 23, 24, 25, 26, 27, 28, 31, 32, 33, 34, 35, 36, 37, 38, 66, 67, 68, 77

library(vegan)
# Ignore snow data, limited history and unreliable
clim_prelim <- full_clim[c(1, 9, 10, 11, 12, 15, 19, 23, 24, 25, 26, 27, 28, 31, 32, 33, 34, 35, 36, 37, 38)]

for (var in colnames(clim_prelim)){
  clim_var <- full_clim[[var]]
  try(prelim_gam <- gam(forcing_QFA~s(clim_var)))
  plot(prelim_gam, xlab=var)
}

# Precipitation data clusters
biplot(princomp(prelim_clim))
arm::corrplot(prelim_clim)
plot(hclust(dist(t(prelim_clim))))

# Relevant subset
# iWUE // more of a mediating variable?
# CO2 // is this better than WUE? spurious correlation?
# Annual minimum temperature
# Mean temperature of period 3
# Total precipitation for period 3
# Number.of.days.of.growing.season

# GAM models on hypothesized relations

# Compiling relevant independent variables
sel_clim <- full_clim[c("Annual.minimum.temperature", "gdd.above.base_temp.for.period.3", "Total.precipitation.for.period.3","Number.of.days.of.growing.season")]
names(sel_clim) <- c("annual_min", "growing_season_GDD", "growing_season_prec", "growing_season_length")
sel_clim$WUE <- F_Wi[as.numeric(names(F_Wi))<=2007 & as.numeric(names(F_Wi))>=1901]
sel_clim$CO2 <- Ca[as.numeric(rownames(Ca))<=2007 & as.numeric(rownames(Ca))>=1901,]
sel_clim$sWUE <- sel_clim$WUE/sel_clim$CO2

# Experiment with simple detrending
# library(pracma)
# 
# forcing_QFA <- pracma::detrend(forcing_QFA)[,1]
# forcing_FA <- pracma::detrend(forcing_FA)[,1]


# Run the models
clim_model_QFA <- gam(log(forcing_QFA)~
      s(sel_clim[["WUE"]])+
#       s(sel_clim[["CO2"]])+
#       s(forcing_QFA_l1)+               
#       s(sel_clim[["annual_min"]])+
      s(sel_clim[["growing_season_GDD"]])+
      s(sel_clim[["growing_season_prec"]])+
      s(sel_clim[["growing_season_length"]]))

clim_model_FA <- gam(log(forcing_FA)~
      s(sel_clim[["WUE"]])+
#       s(sel_clim[["CO2"]])+
#       s(forcing_QFA_l1)+               
#       s(sel_clim[["annual_min"]])+
      s(sel_clim[["growing_season_GDD"]])+
      s(sel_clim[["growing_season_prec"]])+
      s(sel_clim[["growing_season_length"]]))

clim_model_lin_QFA <- lm(log(forcing_QFA)~
     (sel_clim[["WUE"]])+
#      (sel_clim[["CO2"]])+
#      (forcing_QFA_l1)+               
#      (sel_clim[["annual_min"]])+
     (sel_clim[["growing_season_GDD"]])+
     (sel_clim[["growing_season_prec"]])+
     (sel_clim[["growing_season_length"]]))

clim_model_lin_FA <- gam(log(forcing_FA)~
      (sel_clim[["WUE"]])+
#       (sel_clim[["CO2"]])+
#       (forcing_QFA_l1)+               
#       (sel_clim[["annual_min"]])+
      (sel_clim[["growing_season_GDD"]])+
      (sel_clim[["growing_season_prec"]])+
      (sel_clim[["growing_season_length"]]))

summary(clim_model_QFA)
summary(clim_model_FA)

summary (clim_model_lin_QFA)
summary (clim_model_lin_FA)


plot(clim_model_QFA)
plot(clim_model_FA)

# # Modelling WUE and CO2 effects on growth
# co2_growth_QFA <- gam(forcing_QFA~s(sel_clim[["CO2"]]))
# co2_growth_FA <- gam(forcing_FA~s(sel_clim[["CO2"]]))
# 
# wue_growth_QFA <- gam(forcing_QFA~s(sel_clim[["WUE"]]))
# wue_growth_FA <- gam(forcing_FA~s(sel_clim[["WUE"]]))
# 
# # Lagged dependent variables
# forcing_QFA_l1 <- c(forcing_QFA[-1], NA)
# forcing_FA_l1 <- c(forcing_FA[-1], NA)
# 
# lag_QFA <- gam(forcing_QFA~s(forcing_QFA_l1))
# lag_FA <- gam(forcing_FA~s(forcing_FA_l1))
# 
# # Report results of experimentation
# summary(co2_growth_QFA)$r.sq
# summary(co2_growth_FA)$r.sq
# 
# summary(wue_growth_QFA)$r.sq
# summary(wue_growth_FA)$r.sq
# 
# summary(lag_QFA)$r.sq
# summary(lag_FA)$r.sq

# Plotting ####

predictors <- c("WUE", "growing_season_GDD", "growing_season_prec", "growing_season_length")
predictor_names <- c("iWUE", "Growing season GDD", "Growing season precipitation (mm)", "Growing season length (days)")
models <- list(clim_model_QFA, clim_model_FA, clim_model_lin_QFA, clim_model_lin_FA)
model_names <- c("GAM_QFA", "GAM_FA", "LM_QFA", "LM_FA")


partial_response <- data.frame("Model"=NA, "var"=NA, "x"=NA, "response"=NA, "se"=NA)[0,]
for (i in 1:length(models)){
  model <- models[[i]]
  
  fit_se <- predict (model, type="terms", se.fit=T)
  fit <- fit_se$fit
  se <- fit_se$se
  
  colnames(fit) <- predictors
  colnames(se) <- predictors
  
  partial_response_i <- data.frame(melt(fit), melt(se)$value) 
  names(partial_response_i) <- c("x", "var", "response", "se")
  
  partial_response_i$x <- melt(sel_clim[predictors])$value
  
  partial_response_i$Model <- model_names[i]
  
  partial_response <- rbind(partial_response, partial_response_i)  
}

for (v in 1:length(predictors)){
  df <- partial_response[partial_response$var==predictors[v],]
  print(ggplot (data=df, aes(x=x, y=response, colour=Model, fill=Model, ymax = response + 1.96*se, ymin=response - 1.96*se))+geom_line(size=1.5)+labs(y="Partial predictions", x=predictor_names[v])+theme_bw()+geom_ribbon(alpha=0.1)+geom_ribbon(alpha=0.1)+geom_hline(y=0))
}

# Plot the predictor variables vs. time
for (v in predictors){
  df <- data.frame(y=sel_clim[[v]], x=as.numeric(rownames(sel_clim)))
  limits <- aes()
  print(ggplot(data=df, aes(x=x, y=y)) + geom_line(size=1.5)+theme_bw() + xlab("Time") + ylab(predictor_names[which(predictors==v)]))
}

# Exploring WUE vs. climate links ####

wue_clim_fit <- vector()
for (var in 1:ncol(full_clim)){
  clim_var <- full_clim[[var]]
  try(prelim_gam <- gam(sel_clim$sWUE~s(clim_var)))
  wue_clim_fit <- c(wue_clim_fit, summary(prelim_gam)$r.sq)
}
names(wue_clim_fit) <- names(full_clim)
wue_clim_fit <- sort(wue_clim_fit, decreasing=T)
print(wue_clim_fit)

# Scaled WUE seems to be a mediating variable; indicates drought stress, declining as precipitation increasees

clim_model_QFA <- gam(log(forcing_QFA)~
                        s(sel_clim[["WUE"]])+
                        #       s(sel_clim[["CO2"]])+
                        #       s(forcing_QFA_l1)+               
                        #       s(sel_clim[["annual_min"]])+
                        s(sel_clim[["growing_season_GDD"]])+
                        s(sel_clim[["growing_season_prec"]])+
                        s(sel_clim[["growing_season_length"]]))


