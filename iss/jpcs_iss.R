# Load dependencies ####
library(reshape2)
library(dplR)
library(ggplot2)
library(plyr)
library(mgcv)

# Saving options ####

# Root save path
main_wd <- getwd()
figure_path <- paste(getwd(), "Figures", "Pinus_banksiana", sep="/")

# General pdf saving parameters
pdf_width=11
pdf_height=8.5

# Load data ####

jp <- read.rwl ("./Data/jpcs.rwl")

# Remove trees less than 30 years old
# young_trees <- apply(jp, FUN=function(x){sum(!is.na(x))}, MARGIN=2)<=30
# 
# jp <- jp[!young_trees]

# Remove recent years
jp <- jp[-which((rownames(jp) %in% c("2011"))),]

# Standardize chronology ####

# Standardize using modified negative exponential curve
jp_exp <- detrend(jp, method="ModNegExp", make.plot=F)
jp_exp_chron <- chron(jp_exp, prefix="jpx", biweight=F, prewhiten=F)

# Standardize using splines
jp_spl <- detrend(jp, method="Spline", make.plot=F)
jp_spl_chron <- chron(jp_spl, prefix="jps", biweight=F, prewhiten=F)

# Standardize using flat detrending
jp_mean <- detrend(jp, method="Mean", make.plot=F)
jp_mean_chron <- chron(jp_mean, prefix="jpm", biweight=F, prewhiten=F)

# Combining information
jp_chron <- data.frame(Year=as.numeric(rownames(jp_exp_chron)), Exp=jp_exp_chron[[1]], Spl=jp_spl_chron[[1]], Flat =jp_mean_chron[[1]], Samples=jp_spl_chron[[2]])

jp_chron_melt <- melt(jp_chron[,1:4], id.vars="Year")

jp_chron_all_plot <- ggplot(jp_chron_melt, aes(x=Year, y=value, colour=variable)) + geom_line(size=1) + theme_bw() + geom_hline(y=1) + ylab("Standardized chronology value")

jp_chron_exp_plot <- ggplot(jp_chron_melt[jp_chron_melt$variable=="Exp",], aes(x=Year, y=value)) + geom_line() + theme_bw() + geom_hline(y=1) + ylab("Standardized chronology value")

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

# Convert format to a tree ring table
deltaC_Wood <-  dcast(raw_deltaC_Wood, Year~Site, value.var="X.13C.Wood", mean)
deltaC_Wood[sapply(deltaC_Wood, is.nan)] <- NA
rownames(deltaC_Wood) <- deltaC_Wood$Year
deltaC_Wood <- deltaC_Wood[-1]

# Convert to DeltaC
DeltaC <- (deltaC_Atm[[1]]-deltaC_Wood)/(1000+deltaC_Wood)*1000

# Convert to iWUE
# Do not Standardize iWUE data

a <- 4.4
b <- 27
Wi <- Ca[[2]]/1.6*(1-(DeltaC-a)/(b-a))

# Find average iWUE
Wi_chron <- apply(Wi, FUN=mean, MARGIN=1, na.rm=T)

# Plots of delta, Delta and iWUE

d13C_df <- data.frame(Year=as.numeric(rownames(deltaC_Wood)), value=rowMeans(deltaC_Wood, na.rm=T), reliable=(apply(deltaC_Wood, FUN=function(x){sum(!is.na(x))}, MARGIN=1)>=5))

D13C_df <- data.frame(Year=as.numeric(rownames(DeltaC)), value=rowMeans(DeltaC, na.rm=T), reliable=(apply(DeltaC, FUN=function(x){sum(!is.na(x))}, MARGIN=1)>=5))

Wi_df <- data.frame(Year=as.numeric(rownames(DeltaC)), value=Wi_chron, reliable=(apply(Wi, FUN=function(x){sum(!is.na(x))}, MARGIN=1)>=5))


d13C_plot <- ggplot(d13C_df, aes(x=Year, y=value, linetype=reliable)) + geom_line() + scale_x_continuous(breaks=c(1850, 1900, 1950, 2000)) + theme_bw() + ylab ("delta13C") + scale_linetype_manual(values=c("dotted", "solid")) + theme(legend.position="none")

D13C_plot <- ggplot(D13C_df, aes(x=Year, y=value, linetype=reliable)) + geom_line() + scale_x_continuous(breaks=c(1850, 1900, 1950, 2000)) + theme_bw() + ylab ("Delta13C") + scale_linetype_manual(values=c("dotted", "solid")) + theme(legend.position="none")

Wi_plot <- ggplot(Wi_df, aes(x=Year, y=value, linetype=reliable)) + geom_line() + scale_x_continuous(breaks=c(1850, 1900, 1950, 2000)) + theme_bw() + ylab ("iWUE") + scale_linetype_manual(values=c("dotted", "solid")) + theme(legend.position="none")

# Find climate-growth links ####

# Load climate data
jp.clim <- read.csv("./Data/IgnaceClimate.csv")
foo <- melt(jp.clim, id.vars=c("Month", "Year"))
jp.clim.wide <- dcast(foo, Year~variable+Month) 
full_clim <- as.data.frame(read.csv("./Data/full_clim.csv", header=T))
rownames(full_clim) <- full_clim[[1]]

# Remove metadata info
full_clim <- full_clim[-1:-5]

# GAM and generalized linear models for climate response

# Summer drought: June max temp
# Winter prec.: February + March prec
# Summer precipitation
# Growing season length
# Growing season GDD
# WUE

# Compiling relevant independent variables
sel_clim <- full_clim[c("gdd.above.base_temp.for.period.3", "Total.precipitation.for.period.3","Number.of.days.of.growing.season", "June.mean.monthly.maximum.temperature")]
names(sel_clim) <- c("GDD", "summer_prec", "growing_season_length", "june_max")
#sel_clim$winter_prec <- full_clim$February.mean.monthly.precipitation + full_clim$March.mean.monthly.precipitation
sel_clim$year <- as.numeric(rownames(sel_clim))

# Truncate chronology to only matching years
# Pick which detrended chronology to use
rw_response <- jp_chron[which(jp_chron$Year %in% sel_clim$year),"Exp"]

# Truncate D13C chronology to include only reliable years
D13C_chron <- apply(DeltaC, FUN=mean, MARGIN=1, na.rm=T)
D13C_response <- D13C_chron[which(as.numeric(names(D13C_chron))>=1915 & as.numeric(names(D13C_chron))<=2007)]

# Truncate data to exclude unreliable WUE data
D13C_clim <- sel_clim[sel_clim$year >= 1915,]

# Ring width x climate models
clim_rw_glm <- gam(rw_response~
                        sel_clim[["GDD"]]+
                        sel_clim[["june_max"]]+                                    
                        sel_clim[["summer_prec"]]+
                        #sel_clim[["winter_prec"]]+                  
                        sel_clim[["growing_season_length"]], 
                        family=gaussian(link="log"))

clim_rw_gam <- gam(rw_response~
                        s(sel_clim[["GDD"]])+
                        s(sel_clim[["june_max"]])+                                    
                        s(sel_clim[["summer_prec"]])+
                        #s(sel_clim[["winter_prec"]])+                  
                        s(sel_clim[["growing_season_length"]]), 
                        family=gaussian(link="log"))

summary(clim_rw_glm)
summary(clim_rw_gam)


# iWUE x climate models
clim_D13C_glm <- gam(D13C_response~
                     D13C_clim[["GDD"]]+
                     D13C_clim[["june_max"]]+                                    
                     D13C_clim[["summer_prec"]]+
                     #D13C_clim[["winter_prec"]]+                  
                     D13C_clim[["growing_season_length"]], 
                   family=gaussian())

clim_D13C_gam <- gam(D13C_response~
                     s(D13C_clim[["GDD"]])+
                     s(D13C_clim[["june_max"]])+                                    
                     s(D13C_clim[["summer_prec"]])+
                     #s(D13C_clim[["winter_prec"]])+                  
                     s(D13C_clim[["growing_season_length"]]), 
                   family=gaussian())

summary(clim_D13C_glm)
summary(clim_D13C_gam)

# Plotting climate response analysis ####

# Model terms and settings
predictors <- c("GDD", 
                "june_max", 
                "summer_prec", 
                #"winter_prec", 
                "growing_season_length")

predictor_names <- c("Growing season GDD", 
                   "June maximum temperature",  
                   "Summer precipitation (mm)", 
                   #"Winter precipitation (mm)", 
                   "Growing season length (days)")

models <- list(clim_rw_glm, clim_D13C_glm)
model_names <- c("Ring Width", "D13C")

model_years <- list(1901:2007, 1915:2007)

# Finding partial responses for each model
partial_response <- data.frame("Model"=NA, "var"=NA, "x"=NA, "response"=NA, "se"=NA)[0,]
for (i in 1:length(models)){
  model_i <- models[[i]]
  
  fit_se <- predict (model_i, type="terms", se.fit=T)
  fit <- fit_se$fit
  se <- fit_se$se
  
  colnames(fit) <- predictors
  colnames(se) <- predictors
  
  partial_response_i <- data.frame(melt(fit), melt(se)$value) 
  names(partial_response_i) <- c("x", "var", "response", "se")
    
  partial_response_i$data <- melt(sel_clim[which(sel_clim$year %in% model_years[[i]]),][predictors])$value
  
  partial_response_i$year <- model_years[[i]]
  
  partial_response_i$Model <- model_names[i]
  
  partial_response <- rbind(partial_response, partial_response_i)  
}

# Renaming levels for pretty graph label
levels(partial_response$var) <- c("Growing season GDD", 
                                  "Growing season length (days)", 
                                  "June maximum temperature", 
                                  "Summer precipitation (mm)", 
                                  #"Winter precipitation (mm)"
                                  )
# Rearranging levels to control order on plots
# c(5,4,3,1,2,6)
#partial_response$var <- factor(partial_response$var, levels(partial_response$var)[c(4,3,2,1,5)])

# Plot partial predictions

clim_partial_predictions_plot <- ggplot (data=partial_response, aes(x=data, y=response, ymax = response + 1.96*se, ymin=response - 1.96*se))+geom_line(size=1.5)+labs(y="Standardized partial predictions", x="Predictor variable")+theme_bw()+geom_ribbon(alpha=0.3)+geom_hline(y=0) + facet_grid(Model~var, scales="free_x") + theme_bw()

# Plot the impact of each predictor variable by time
clim_predictions_by_time_plot <- ggplot(data=partial_response, aes(x=year, y=response, ymax = response + 1.96*se, ymin=response - 1.96*se)) + facet_grid(Model~var) + geom_line(size=1) + geom_ribbon(alpha=0.3) + geom_hline(y=0) + ylim(c(-1.85, 1.85)) + theme_bw() + labs(y="Standardized partial predictions", x="Year")


# Plot the predictor variables vs. time
clim_by_time_plot <- ggplot(data=partial_response, aes(x=year, y=data)) + facet_grid(var~., scales="free_y") + geom_line(size=1) + theme_bw() + labs(y="Climate variable", x="Year")

# Saving results ####

# Standardization
ggsave("./iss/Figures/jp_chron_all_plot.svg", jp_chron_all_plot, width=8.5, height=3)
ggsave("./iss/Figures/jp_chron_all_plot.pdf", jp_chron_all_plot, width=8.5, height=3)

ggsave("./iss/Figures/jp_chron_exp_plot.svg", jp_chron_exp_plot, width=8.5, height=3)
ggsave("./iss/Figures/jp_chron_exp_plot.pdf", jp_chron_exp_plot, width=8.5, height=3)

# Delta 13C
ggsave("./iss/Figures/d13C_plot.svg", d13C_plot, width=8.5, height=3)
ggsave("./iss/Figures/d13C_plot.pdf", d13C_plot, width=8.5, height=3)

ggsave("./iss/Figures/Del13C_plot.svg", D13C_plot, width=8.5, height=3)
ggsave("./iss/Figures/Del13C_plot.pdf", D13C_plot, width=8.5, height=3)

ggsave("./iss/Figures/Wi_plot.svg", Wi_plot, width=8.5, height=3)
ggsave("./iss/Figures/Wi_plot.pdf", Wi_plot, width=8.5, height=3)

# Climate response
ggsave("./iss/Figures/clim_partial_predictions_plot .svg", clim_partial_predictions_plot, width=8.5, height=3)
ggsave("./iss/Figures/clim_partial_predictions_plot.pdf", clim_partial_predictions_plot, width=8.5, height=3)

ggsave("./iss/Figures/clim_predictions_by_time_plot.svg", clim_predictions_by_time_plot, width=8.5, height=3)
ggsave("./iss/Figures/clim_predictions_by_time_plot.pdf", clim_predictions_by_time_plot, width=8.5, height=3)

ggsave("./iss/Figures/clim_by_time_plot.svg", clim_by_time_plot, width=3, height=8.5)
ggsave("./iss/Figures/clim_by_time_plot.pdf", clim_by_time_plot, width=3, height=8.5)