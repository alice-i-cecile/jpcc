# Load dependencies ####

library(pastecs)
library(fBasics)
library(bbmle)
library(reshape2)
library(dplR)
library(ggplot2)
library(gtools)
library(plyr)

# Saving options ####

# Root save path
main_wd <- getwd()
figure_path <- paste(getwd(), "Figures", "Pinus_banksiana", sep="/")

# General pdf saving parameters
pdf_width=11
pdf_height=8.5

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
  I.names <- dimnames(tra)[[1]]
  
  # Get site code for each core
  sites <- getSite(I.names)
  
  # Find first year of each core
  firstYears <- vector ()
  for (i in I.names){
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

# Load data ####

jp <- read.rwl ("./Data/jpcs.rwl")

# Convert 0 values to NA
jp[which(jp==0, arr.ind=T)] <- NA

# Standardize using modified negative exponential curve

# Standardize using splines

# Plotting sample depth ####
sample_depth_melt <- melt(jp.depth.T)
sample_depth_melt$Year <- as.numeric(rownames(sample_depth_melt))

sample_depth_plot <- ggplot(sample_depth_melt, aes(x=Year, y=value)) + geom_area() + theme_bw() + ylab("Number of series in chronology")

# Plotting effects ####

# Plot I
# Plot I vs age
I_vs_age_plot <- ggplot(I_df, aes(x=age, y=value)) + geom_point() + geom_smooth(colour="black") + facet_grid(variable~.) + theme_bw() + ylab("Individual effect (I)") + xlab("Age")

# Plot histogram of I
I_hist <- ggplot(I_df, aes(x=value, colour=variable, fill=variable)) + geom_density(alpha=0.5)+theme_bw()+xlab("Individual effect (I)") + ylab("Density") + scale_fill_manual(values=c("grey20", "grey50")) + scale_colour_manual(values=c("grey20", "grey50"))+theme(legend.title=element_blank()) + theme(legend.position=c(0.9,0.85)) +  scale_x_log10(breaks=c(0.25, 0.5, 1, 2, 4), lim=c(0.25, 4))+geom_vline(x=1)

# Plot T
T_plot <- ggplot(T_df, aes(x=as.numeric(t), y=value, linetype=reliable)) + geom_line()+ facet_grid(variable~.)+theme_bw()+xlab("Year") + ylab("Time effect (T)") +geom_hline(y=1) + scale_x_continuous(breaks=c(1850, 1900, 1950, 2000)) + scale_linetype_manual(values=c("dotted", "solid")) + theme(legend.position="none")

# Plot A
A_plot <- ggplot(A_df, aes(x=as.numeric(a), y=value, linetype=reliable)) + geom_line()+ facet_grid(variable~.)+theme_bw()+xlab("Age") + ylab("Age effect (A) in mm") + scale_linetype_manual(values=c("dotted", "solid")) + theme(legend.position="none")

# Carbon-13 analysis ####
raw_deltaC_Wood <- read.csv("./Data/tree_13c.csv")

Ca <- read.csv("./Data/co2_record.csv")
Ca <- Ca[Ca$Year >= 1877 & Ca$Year <= 2010, -3]
rownames(Ca) <- Ca$Year

deltaC_Atm <- read.csv("./Data/13C_record.csv")
rownames(deltaC_Atm) <- deltaC_Atm[[1]]
deltaC_Atm <- deltaC_Atm[-1]
deltaC_Atm <- data.frame(DC=deltaC_Atm[as.numeric(rownames(deltaC_Atm))<=2010 & as.numeric(rownames(deltaC_Atm))>=1877,])
rownames(deltaC_Atm) <- 1877:2010

# # Convert format to a tree ring table
deltaC_Wood <-  dcast(raw_deltaC_Wood, Year~Site, value.var="X.13C.Wood", mean)
deltaC_Wood[sapply(deltaC_Wood, is.nan)] <- NA
rownames(deltaC_Wood) <- deltaC_Wood$Year
deltaC_Wood <- deltaC_Wood[-1]

# Convert to DeltaC
DeltaC <- (deltaC_Atm[[1]]-deltaC_Wood)/(1000+deltaC_Wood)*1000

# Convert to iWUE
# Do not Standardize iWUE data

b <- 27
Wi <- Ca[[2]]/1.6*(1-(DeltaC-a)/(b-a))


# Plots of delta and Delta

d13C_plot <- ggplot(data.frame(Year=as.numeric(rownames(deltaC_Wood)), value=rowMeans(deltaC_Wood, na.rm=T), reliable=Wi.depth.T >= 5), aes(x=Year, y=value, linetype=reliable)) + geom_line() + scale_x_continuous(breaks=c(1850, 1900, 1950, 2000)) + theme_bw() + ylab ("delta13C") + scale_linetype_manual(values=c("dotted", "solid")) + theme(legend.position="none")

D13C_plot <- ggplot(data.frame(Year=as.numeric(rownames(DeltaC)), value=rowMeans(DeltaC, na.rm=T), reliable=Wi.depth.T >= 5), aes(x=Year, y=value, linetype=reliable)) + geom_line() + scale_x_continuous(breaks=c(1850, 1900, 1950, 2000)) + theme_bw() + ylab ("Delta13C") + scale_linetype_manual(values=c("dotted", "solid")) + theme(legend.position="none")

# Find climate-growth links ####
# Preprocess forcing info from tree rings
time_ITA <- T_jp$ITA[as.numeric(rownames(T_jp))<=2007 & as.numeric(rownames(T_jp))>=1901]
time_IT <- T_jp$IT[as.numeric(rownames(T_jp))<=2007 & as.numeric(rownames(T_jp))>=1901]
time_TA <- T_jp$TA[as.numeric(rownames(T_jp))<=2007 & as.numeric(rownames(T_jp))>=1901]
time_T <- T_jp$T[as.numeric(rownames(T_jp))<=2007 & as.numeric(rownames(T_jp))>=1901]

names(time_ITA) <- 1901:2007
names(time_IT) <- 1901:2007
names(time_TA) <- 1901:2007
names(time_T) <- 1901:2007

# Load climate data
jp.clim <- read.csv("./Data/IgnaceClimate.csv")
foo <- melt(jp.clim, id.vars=c("Month", "Year"))
jp.clim.wide <- dcast(foo, Year~variable+Month) 
full_clim <- as.data.frame(read.csv("./Data/full_clim.csv", header=T))
rownames(full_clim) <- full_clim[[1]]
full_clim <- full_clim[-1:-5]

# Basic climate correlations and response ####
library(bootRes)

# Computing responses
resp_ITA <- dcc(data.frame(chron=time_ITA, samp.depth=jp.depth.T[30:136]), jp.clim, method="r", start=2, end=10)
resp_IT <- dcc(data.frame(chron=time_IT, samp.depth=jp.depth.T[30:136]), jp.clim, method="r", start=2, end=10)
resp_TA <- dcc(data.frame(chron=time_TA, samp.depth=jp.depth.T[30:136]), jp.clim, method="r", start=2, end=10)
resp_T <- dcc(data.frame(chron=time_T, samp.depth=jp.depth.T[30:136]), jp.clim, method="r", start=2, end=10)

# Cleaning responses
annotate_response <- function(dcc_out, id)
{
  split_code <- strsplit(rownames(dcc_out), split=".curr.")
  dcc_out$var <- sapply(split_code, function(x){x[1]})
  dcc_out$month <- sapply(split_code, function(x){x[2]})
  dcc_out$id <- id
  
  return(dcc_out)
}

resp_ITA <- annotate_response(resp_ITA, "ITA")
resp_IT <- annotate_response(resp_IT, "IT")
resp_TA <- annotate_response(resp_TA, "TA")
resp_T <- annotate_response(resp_T, "T")

resp_all <- rbind(resp_ITA, resp_IT, resp_TA, resp_T)
resp_all$month[resp_all$month=="apr"] <- "April"
resp_all$month[resp_all$month=="aug"] <- "August"
resp_all$month[resp_all$month=="dec"] <- "December"
resp_all$month[resp_all$month=="feb"] <- "February"
resp_all$month[resp_all$month=="jan"] <- "January"
resp_all$month[resp_all$month=="jul"] <- "July"
resp_all$month[resp_all$month=="jun"] <- "June"
resp_all$month[resp_all$month=="mar"] <- "March"
resp_all$month[resp_all$month=="may"] <- "May"
resp_all$month[resp_all$month=="nov"] <- "November"
resp_all$month[resp_all$month=="oct"] <- "October"
resp_all$month[resp_all$month=="sep"] <- "September"

resp_all$month <- factor(resp_all$month, unique(resp_all$month))
resp_all$id <- factor(resp_all$id, unique(resp_all$id))
resp_all$var <- factor(resp_all$var, unique(resp_all$var))

# Plotting
response_plot <- ggplot(resp_all, aes(x=id, y=coef, fill=id)) + geom_bar(stat="identity") + facet_grid(month~var) + theme_bw() + xlab("Standardization model") + ylab("Standardized growth-climate response") + geom_hline(y=0) + geom_linerange(aes(ymin=ci.lower, ymax=ci.upper)) + scale_fill_manual(values=c("darkorchid3", "red", "blue", "grey40")) + guides(fill=FALSE)

# GAM and linear models on prelim relations ####

# Summer drought: June max temp
# Winter prec.: February + March prec
# Summer prec
# Growing season length
# Growing season GDD
# CO2
# WUE

# Compiling relevant independent variables
sel_clim <- full_clim[c("gdd.above.base_temp.for.period.3", "Total.precipitation.for.period.3","Number.of.days.of.growing.season", "June.mean.monthly.maximum.temperature")]
names(sel_clim) <- c("GDD", "summer_prec", "growing_season_length", "june_max")
sel_clim$winter_prec <- full_clim$February.mean.monthly.precipitation + full_clim$March.mean.monthly.precipitation
sel_clim$WUE <- T_Wi[as.numeric(names(T_Wi))<=2007 & as.numeric(names(T_Wi))>=1901]
sel_clim$CO2 <- Ca[as.numeric(rownames(Ca))<=2007 & as.numeric(rownames(Ca))>=1901, 2]
sel_clim$DeltaC <- rowMeans(DeltaC, na.rm=T)[which(names(rowMeans(DeltaC, na.rm=T))==1901):which(names(rowMeans(DeltaC, na.rm=T))==200)]
sel_clim$year <- as.numeric(rownames(sel_clim))

# Truncate data to exclude unreliable WUE data
sel_clim <- sel_clim[sel_clim$year >= 1915,]
sel_T_ITA <- time_ITA[as.numeric(names(time_ITA)) >= 1915]
sel_T_IT <- time_IT[as.numeric(names(time_IT)) >= 1915]
sel_T_TA <- time_TA[as.numeric(names(time_TA)) >= 1915]
sel_T_T <- time_T[as.numeric(names(time_T)) >= 1915]

# Run the models
clim_model_ITA <- gam(log(sel_T_ITA)~
                        s(sel_clim[["WUE"]])+
                        #    s(sel_clim[["CO2"]])+
                        s(sel_clim[["GDD"]])+
                        s(sel_clim[["june_max"]])+                                    
                        #    s(sel_clim[["summer_prec"]])+
                        s(sel_clim[["winter_prec"]])+                  
                        s(sel_clim[["growing_season_length"]]))

clim_model_IT <- gam(log(sel_T_IT)~
                       s(sel_clim[["WUE"]])+
                       #    s(sel_clim[["CO2"]])+
                       s(sel_clim[["GDD"]])+
                       s(sel_clim[["june_max"]])+                                    
                       #    s(sel_clim[["summer_prec"]])+
                       s(sel_clim[["winter_prec"]])+                  
                       s(sel_clim[["growing_season_length"]]))

clim_model_TA <- gam(log(sel_T_TA)~
                       s(sel_clim[["WUE"]])+
                       #    s(sel_clim[["CO2"]])+
                       s(sel_clim[["GDD"]])+
                       s(sel_clim[["june_max"]])+                                    
                       #    s(sel_clim[["summer_prec"]])+
                       s(sel_clim[["winter_prec"]])+                  
                       s(sel_clim[["growing_season_length"]]))

clim_model_T <- gam(log(sel_T_T)~
                      s(sel_clim[["WUE"]])+
                      #     s(sel_clim[["CO2"]])+
                      s(sel_clim[["GDD"]])+
                      s(sel_clim[["june_max"]])+                                    
                      #     s(sel_clim[["summer_prec"]])+
                      s(sel_clim[["winter_prec"]])+                  
                      s(sel_clim[["growing_season_length"]]))

clim_model_lin_ITA <- lm(log(sel_T_ITA)~
                           (sel_clim[["WUE"]])+
                           #     (sel_clim[["CO2"]])+
                           (sel_clim[["GDD"]])+
                           (sel_clim[["june_max"]])+                                    
                           #     (sel_clim[["summer_prec"]])+
                           (sel_clim[["winter_prec"]])+                  
                           (sel_clim[["growing_season_length"]]))

clim_model_lin_IT <- lm(log(sel_T_IT)~
                          (sel_clim[["WUE"]])+
                          #     (sel_clim[["CO2"]])+
                          (sel_clim[["GDD"]])+
                          (sel_clim[["june_max"]])+                                    
                          #     (sel_clim[["summer_prec"]])+
                          (sel_clim[["winter_prec"]])+                  
                          (sel_clim[["growing_season_length"]]))

clim_model_lin_TA <- lm(log(sel_T_TA)~
                          (sel_clim[["WUE"]])+
                          #     (sel_clim[["CO2"]])+
                          (sel_clim[["GDD"]])+
                          (sel_clim[["june_max"]])+                                    
                          #     (sel_clim[["summer_prec"]])+
                          (sel_clim[["winter_prec"]])+                  
                          (sel_clim[["growing_season_length"]]))

clim_model_lin_T  <- lm(log(sel_T_T)~
                          (sel_clim[["WUE"]])+
                          #     (sel_clim[["CO2"]])+
                          (sel_clim[["GDD"]])+
                          (sel_clim[["june_max"]])+                                    
                          #     (sel_clim[["summer_prec"]])+
                          (sel_clim[["winter_prec"]])+                  
                          (sel_clim[["growing_season_length"]]))

summary(clim_model_ITA)
summary(clim_model_IT)
summary(clim_model_TA)
summary(clim_model_T)

summary (clim_model_lin_ITA)
summary (clim_model_lin_IT)
summary (clim_model_lin_TA)
summary (clim_model_lin_T)

# plot(clim_model_ITA)
# plot(clim_model_IT)
# plot(clim_model_TA)
# plot(clim_model_T)

Rsq_df <- data.frame(
  ITA=c(summary(clim_model_ITA)$r.sq, summary(clim_model_lin_ITA)$r.sq),
  IT=c(summary(clim_model_IT)$r.sq, summary(clim_model_lin_IT)$r.sq),
  TA=c(summary(clim_model_TA)$r.sq, summary(clim_model_lin_TA)$r.sq),
  T=c(summary(clim_model_T)$r.sq, summary(clim_model_lin_T)$r.sq) 
)

rownames(Rsq_df) <- c("gam", "lm")

AIC_df <- data.frame(AIC=c(AIC(clim_model_ITA), AIC(clim_model_lin_ITA)), BIC=c(BIC(clim_model_ITA), BIC(clim_model_lin_ITA)))
rownames(AIC_df) <- c("gam", "lm")

AIC_df$AIC <- AIC_df$AIC - AIC_df$AIC[1]
AIC_df$BIC <- AIC_df$BIC - AIC_df$BIC[1]

# Plotting climate response analysis ####

predictors <- c("WUE", "GDD", "june_max", 
                #"summer_prec", 
                "winter_prec", "growing_season_length")
predictor_names <- c("iWUE", "Growing season GDD", "June maximum temperature",  
                     #"Summer precipitation (mm)", 
                     "Winter precipitation (mm)", "Growing season length (days)")

models <- list(clim_model_ITA, clim_model_lin_ITA)
model_names <- c("GAM", "LM")

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

levels(partial_response$var) <- c("Growing season GDD", 
                                  "Growing season length (days)", 
                                  "June maximum temperature", 
                                  #"Summer precipitation (mm)", 
                                  "Winter precipitation (mm)", 
                                  "iWUE")
# c(5,4,3,1,2,6)
partial_response$var <- factor(partial_response$var, levels(partial_response$var)[c(4,3,2,1,5)])

# Plot partial predictions

clim_model_plot <- ggplot (data=partial_response, aes(x=x, y=response, colour=Model, fill=Model, ymax = response + 1.96*se, ymin=response - 1.96*se))+geom_line(size=1.5)+labs(y="Standardized partial predictions", x="Predictor variable")+theme_bw()+geom_ribbon(alpha=0.1)+geom_hline(y=0) + facet_wrap(~var, scales="free_x", ncol=1) + ylim(c(-0.4, 0.4)) + theme(legend.position="top")

# Plot the impact of each predictor variable by time
impact_by_time <- data.frame(Year=NA, variable=NA, value=NA, model=NA)[0,]
for (i in 1:length(models)){
  model <- models[[i]]
  
  partial_predictions <- as.data.frame(predict (model, type="terms"))
  colnames(partial_predictions) <- c("iWUE", 
                                     "Growing season GDD", 
                                     "June maximum temperature",  
                                     #"Summer precipitation (mm)", 
                                     "Winter precipitation (mm)", 
                                     "Growing season length (days)")
  # partial_predictions$intercept <- model$coefficients[1]
  partial_predictions$year <- rownames(partial_predictions)
  
  partial_predictions <- melt(partial_predictions)
  partial_predictions$year <- as.numeric(partial_predictions$year)
  partial_predictions$model <- c("GAM", "LM")[i]
  
  impact_by_time <- rbind(impact_by_time, partial_predictions)
}

# Develop a composite prediction, weighted by AIC
wide_ibt_lm <- dcast(impact_by_time[impact_by_time$model=="LM",], year~variable)
wide_ibt_gam <- dcast(impact_by_time[impact_by_time$model=="GAM",], year~variable)

lm_weight <- exp(-0.5*(AIC(clim_model_lin_ITA)-AIC(clim_model_ITA)))
gam_weight <- 1-lm_weight

wide_ibt_composite <- wide_ibt_lm * lm_weight + wide_ibt_gam * gam_weight
long_ibt_composite <- melt(wide_ibt_composite, id.vars="year")
long_ibt_composite$model <- "Composite"

impact_by_time <- rbind(impact_by_time, long_ibt_composite)

# Plotting
impact_by_time$model <- factor(impact_by_time$model, levels=unique(impact_by_time$model)[c(2,1,3)])
impact_by_time$variable <- factor(impact_by_time$variable, levels=unique(impact_by_time$variable)[c(4,3,5,2,1)])

impact_by_time_plot <- ggplot(impact_by_time, aes(x=year, y=value)) + geom_line() + facet_grid(variable~model) +theme_bw() + ylab("Partial predictions") + xlab("Year")


# Plot the predictor variables vs. time
clim_df <- sel_clim
names(clim_df) <- c("iWUE", "Growing season GDD", "Summer precipitation (mm)", "Winter precipitation (mm)", "Growing season length (days)", "Ambient CO2 (ppm)")
clim_df$Year <- as.numeric(rownames(sel_clim))
clim_df <- melt(clim_df, id.vars="Year")
clim_df$variable <- factor(clim_df$variable, levels(clim_df$variable)[c(4,3,2,5,1,6)])
# [c(5,4,3,2,6,1,7)]

clim_patterns_plot <- ggplot(data=clim_df, aes(x=Year, y=value)) + geom_line(size=1.5)+theme_bw() + xlab("Time") + ylab("Predictor variable") + facet_wrap(~variable, scales="free_y", ncol=1)

# Exploring WUE vs. climate links ####

# First, we need to remove the CO2 signal
# Saurer et al. 2004 says that Ci/Ca should be roughly constant
# Carbon isotope discrimination indicates improving water-use efﬁciency of trees in northern Eurasia over the last 100 years
# DeltaC = a + (b-a)Ci/Ca
# Wi = Ca[[2]]/1.6*(1-(DeltaC-a)/(b-a))
# Wi = Ca/1.6*(1-(a + (b-a)Ci/Ca-a)/(b-a))
# Wi = Ca/1.6*(1-Ci/Ca)
# Wi = Ca/1.6-Ci/1.6
# Wi = (Ca-Ci)/1.6
# Ci = Ca-1.6*Wi

# Wi = (1-B_Ca)*Ca + B_clim*Clim + noise

Wi_Ca_model <- gam(1.6*WUE ~ CO2 + 0, data=sel_clim)

Wi_lin_model <- gam(1.6*WUE 
                    ~ CO2 +
                      (sel_clim[["GDD"]])+
                      (sel_clim[["june_max"]])+                                    
                      #(sel_clim[["summer_prec"]])+
                      (sel_clim[["winter_prec"]])+                  
                      (sel_clim[["growing_season_length"]])+
                      0, 
                    data=sel_clim)

Wi_gam_model <- gam(1.6*WUE 
                    ~ CO2 +
                      s(sel_clim[["GDD"]])+
                      s(sel_clim[["june_max"]])+                                    
                      #s(sel_clim[["summer_prec"]])+
                      s(sel_clim[["winter_prec"]])+                  
                      s(sel_clim[["growing_season_length"]])+
                      0, 
                    data=sel_clim)

summary(Wi_Ca_model)
summary(Wi_lin_model)
summary(Wi_gam_model)

plot(Wi_gam_model)

Wi_Delta_model <- lm(DeltaC 
                     ~ (sel_clim[["GDD"]])+
                       (sel_clim[["june_max"]])+                                    
                       #(sel_clim[["summer_prec"]])+
                       (sel_clim[["winter_prec"]])+                  
                       (sel_clim[["growing_season_length"]]), 
                     data=sel_clim)


# Saving results ####

# Sample depth
ggsave("./Figures/sample_depth.svg", sample_depth_plot, width=8.5, height=3)
ggsave("./Figures/sample_depth.pdf", sample_depth_plot, width=8.5, height=3)

# Standardization
save(jp.sfs.ITA, file="./Results/jp.sfs.ITA.RData")
save(jp.sfs.ITA, file="./Results/jp.sfs.IT.RData")
save(jp.sfs.TA, file="./Results/jp.sfs.TA.RData")
save(jp.sfs.T, file="./Results/jp.sfs.T.RData")

write.csv(jp.fit, file="./Results/jp_fit.csv")

ggsave("./Figures/I_vs_age.svg", I_vs_age_plot, width=8.5, height=3)
ggsave("./Figures/I_vs_age.pdf", I_vs_age_plot, width=8.5, height=3)
ggsave("./Figures/I_hist.svg", I_hist, width=8.5, height=3)
ggsave("./Figures/I_hist.pdf", I_hist, width=8.5, height=3)
ggsave("./Figures/T_plot.svg", T_plot, width=8.5, height=6)
ggsave("./Figures/T_plot.pdf", T_plot, width=8.5, height=6)
ggsave("./Figures/A_plot.svg", A_plot, width=8.5, height=3)
ggsave("./Figures/A_plot.pdf", A_plot, width=8.5, height=3)

save(I_age_corr, file="./Results/I_age_corr.RData")
save(T_trend, file="./Results/T_trend.RData")
save(A_trend, file="./Results/A_trend.RData")

write.csv(cor_I, file="./Results/I_cor.csv")
write.csv(cor_T, file="./Results/T_cor.csv")
write.csv(cor_A, file="./Results/A_cor.csv")

# WUE standardization
write.csv(DeltaC, file="./Results/DeltaC.csv")
write.csv(Wi, file="./Results/iWUE.csv")

save(Wi_fes, file="./Results/Wi_fes.RData")
save(Wi_mean, file="./Results/Wi_mean.RData")

write.csv(Wi_fit, file="./Results/Wi_fit.csv")

d13C_plot
ggsave("./Figures/d13C_plot.svg", d13C_plot, width=8.5, height=3)
ggsave("./Figures/d13C_plot.pdf", d13C_plot, width=8.5, height=3)
ggsave("./Figures/D13C_plot.svg", D13C_plot, width=8.5, height=3)
ggsave("./Figures/D13C_plot.pdf", D13C_plot, width=8.5, height=3)


ggsave("./Figures/Wi_I_hist.svg", Wi_I_hist, width=8.5, height=3)
ggsave("./Figures/Wi_I_hist.pdf", Wi_I_hist, width=8.5, height=3)
ggsave("./Figures/Wi_T_plot.svg", Wi_T_plot, width=8.5, height=6)
ggsave("./Figures/Wi_T_plot.pdf", Wi_T_plot, width=8.5, height=6)
ggsave("./Figures/Wi_A_plot.svg", Wi_A_plot, width=8.5, height=3)
ggsave("./Figures/Wi_A_plot.pdf", Wi_A_plot, width=8.5, height=3)

# Linear bootstrapped response
write.csv(resp_all, file="./Results/bootstrapped_clim_responses.csv")
ggsave("./Figures/response_mega_plot.svg", response_plot, width=8.5, height=11)
ggsave("./Figures/response_mega_plot.pdf",response_plot, width=8.5, height=11)

# Climate patterns
ggsave("./Figures/clim_patterns_plot.svg", clim_patterns_plot, width=5, height=11)
ggsave("./Figures/clim_patterns_plot.pdf", clim_patterns_plot, width=5, height=11)

# Climate response
summ_clim_model_ITA <- summary(clim_model_ITA)
summ_clim_model_lin_ITA <- summary(clim_model_lin_ITA)
save(summ_clim_model_ITA, file="./Results/summ_clim_model_ITA.RData")
save(summ_clim_model_lin_ITA, file="./Results/summ_clim_model_lin_ITA.RData")

write.csv(Rsq_df, file="./Results/Rsq_clim_response_models.csv")
write.csv(AIC_df, file="./Results/AIC_clim_response_models.csv")

ggsave("./Figures/clim_model_plot.svg", clim_model_plot, width=5, height=11)
ggsave("./Figures/clim_model_plot.pdf", clim_model_plot, width=5, height=11)

ggsave("./Figures/impact_by_time_plot.svg", impact_by_time_plot, width=5, height=11)
ggsave("./Figures/impact_by_time_plot.pdf", impact_by_time_plot, width=5, height=11)