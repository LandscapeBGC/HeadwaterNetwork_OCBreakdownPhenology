#### All packages required for Headwater Network OC Breakdown Flux Analysis
#### D.Hare
#### danielle.hare@uconn.edu


library(tidyverse)
library(lubridate)
library(zoo)
library(mgcv)
library(splines2)
library(DVstats) #if errors occur with dependancy packages - open a new session and download there - seems to work... 
library(dataRetrieval)
library(smwrBase)
library(baytrends) #fillMissing for dates
library(gridExtra)
library(grid)
library(kableExtra)
library(lme4)
library(sf)
library(sp)
library(timetk)
library(baytrends) #for fillMissing
library(igraph)
library(broom)
library(sjPlot)
library(performance)
library(lattice) #qqplots
library(cowplot)

SinPara <- function(temp){
  time_units <- 365 #units_day 
  temp$timex <- 2*pi*yday(temp$date)/time_units #convert to radians,yday is yearday
  
  #fit raw data to sin
  temp.lm <- lm(temp$tempD ~ sin(temp$timex) + cos(temp$timex)) #dont know why pipe didnt work...
  
  #equation fit added to dataframe
  temp$Tsin.lm <- coef(temp.lm)['sin(temp$timex)']*sin(temp$timex) +
    coef(temp.lm)['cos(temp$timex)']*cos(temp$timex) +
    coef(temp.lm)['(Intercept)']
  
  #Calculate Phase in Days
  phase <- (time_units/(2*pi))*((3*pi/2) - atan(coef(temp.lm)['cos(temp$timex)']/
                                                  coef(temp.lm)['sin(temp$timex)']))  # units of days
  #Calculate Amplitude
  amp <- sqrt((coef(temp.lm)['sin(temp$timex)']^2) + (coef(temp.lm)['cos(temp$timex)']^2))
  
  #Calculate Mean 
  ymean <- coef(temp.lm)['(Intercept)']
  
  #extract adjusted r2 of regression model
  Rsq_adj <-  round(summary(temp.lm)$adj.r.squared,2)
  #extract residual standard error of regression model
  RSE <- summary(temp.lm)$sigma
  ID <- unique(temp$ID)
  
  #add another WS_ID to make sure it is mainatined (new name)
  sinpara_result <- c(phase, amp, ymean, RSE, Rsq_adj)
  names(sinpara_result) <-  c('phase','amp', 'ymean', "RSE", "Rsq_adj") #simply names
  result <- list(temp, sinpara_result)
  
  return(result)
}


lit_breakdown_pred_BM <-  function(model, input_df){
  require(merTools)
  #uncertainity in the temperature parameter
  #pred_AcSh <- exp(predictInterval(model, newdata = input_df, level = 0.95, which = c("fixed"), include.resid.var=0, ignore.fixed.terms = c("mean_flow_st")))
  
  pred_k <- exp(predict(model, newdata=input_df, re.form = ~0))
  
  b <- bootMer(model, nsim=500, 
               FUN=function(x)predict(x, newdata=input_df, re.form = ~0))
  
  predCL <- exp(t(apply(b$t, MARGIN = 2, FUN = quantile, probs = c(0.025, 0.975))))
  
  return(cbind(pred_k, predCL))
  
}

fragmentation_k_calc <- function(input_df, scaled_TQ_df, mod_k_AF, mod_k_RF){
  #transform tempC and discharge to scaled values 
  Boltz <- (8.617*(10^-5))
  #Acer Data
  input_df_AF <- input_df %>%
    mutate( one.k.T.cent = (1/((tempC + 273.15)* Boltz)) - dplyr::filter(scaled_TQ_df, leaf == "A" & type =="Sh")$avgT,
            mean_flow_st =  if_else(Qout_ls > 1500, scale(1500, center = dplyr::filter(scaled_TQ_df, leaf == "A" & type =="Sh")$avgQ, 
                                                          scale = dplyr::filter(scaled_TQ_df, leaf == "A" & type =="Sh")$stdevQ), 
                                    scale(Qout_ls, center = dplyr::filter(scaled_TQ_df, leaf == "A" & type =="Sh")$avgQ, 
                                          scale = dplyr::filter(scaled_TQ_df, leaf == "A" & type =="Sh")$stdevQ))
    )
  #Rhodo Data
  input_df_RF <- input_df %>%
    mutate( 
      one.k.T.cent = (1/((tempC + 273.15)* Boltz)) - dplyr::filter(scaled_TQ_df, leaf == "R" & type =="Sh")$avgT,
      mean_flow_st =  if_else(Qout_ls > 1500, scale(1500, center = dplyr::filter(scaled_TQ_df, leaf == "R" & type =="Sh")$avgQ, 
                                                    scale = dplyr::filter(scaled_TQ_df, leaf == "R" & type =="Sh")$stdevQ), 
                              scale(Qout_ls, center = dplyr::filter(scaled_TQ_df, leaf == "R" & type =="Sh")$avgQ, 
                                    scale = dplyr::filter(scaled_TQ_df, leaf == "R" & type =="Sh")$stdevQ))
    )
  
  #BM choosen to exclude random effects from dataframe, thus no need for varying deployment days/years and streams 
  df <- distinct(input_df_AF) #to reduce processing time 
  k_ci <- as.data.frame(lit_breakdown_pred_BM(mod_k_AF, df))%>%
    rename("kAF" = pred_k,
           "kAF_UCI" = `97.5%`,
           "kAF_LCI" = `2.5%`)
  
  output_kAF <- cbind(df, k_ci) %>%
    left_join(input_df, ., by = c("date", "tempC", "Qout_ls")) 
  
  #make a dataframe with only unique input parameters for Rhodo Data
  df <- distinct(input_df_RF)
  k_ci <- as.data.frame(lit_breakdown_pred_BM(mod_k_RF, df))
  
  output_kF <- cbind(df, k_ci) %>%
    left_join(output_kAF, ., by = c("date", "tempC", "Qout_ls")) %>% #join with kAF values as well
    rename("kRF" = pred_k,
           "kRF_UCI" = `97.5%`,
           "kRF_LCI" = `2.5%`)
  
  return(output_kF)
  
}


microbial_k_calc <- function(input_df, scaled_TQ_df, mod_k_AM, mod_k_RM){
  Boltz <- (8.617*(10^-5))
  input_df_AM <- input_df %>%
    mutate( one.k.T.cent = (1/((tempC + 273.15)* Boltz)) - dplyr::filter(scaled_TQ_df, leaf == "A" & type =="M")$avgT,
            mean_flow_st =  if_else(Qout_ls > 1500, scale(1500, center = dplyr::filter(scaled_TQ_df, leaf == "A" & type =="M")$avgQ, 
                                                          scale = dplyr::filter(scaled_TQ_df, leaf == "A" & type =="M")$stdevQ), 
                                    scale(Qout_ls, center = dplyr::filter(scaled_TQ_df, leaf == "A" & type =="M")$avgQ, 
                                          scale = dplyr::filter(scaled_TQ_df, leaf == "A" & type =="M")$stdevQ))
    )
  #Rhodo Data
  input_df_RM <- input_df %>%
    mutate( 
      one.k.T.cent = (1/((tempC + 273.15)* Boltz)) - dplyr::filter(scaled_TQ_df, leaf == "R" & type =="M")$avgT,
      mean_flow_st =  if_else(Qout_ls > 1500, scale(1500, center = dplyr::filter(scaled_TQ_df, leaf == "R" & type =="M")$avgQ, 
                                                    scale = dplyr::filter(scaled_TQ_df, leaf == "R" & type =="M")$stdevQ), 
                              scale(Qout_ls, center = dplyr::filter(scaled_TQ_df, leaf == "R" & type =="M")$avgQ, 
                                    scale = dplyr::filter(scaled_TQ_df, leaf == "R" & type =="M")$stdevQ))
    )
  
  
  df_AM <- distinct(input_df_AM) #to reduce processing time 
  df_AM$kAM <- exp(predict(mod_k_AM, newdata=df_AM, re.form = ~0))
  
  #make a dataframe with only unique input parameters for Rhodo Data
  df_RM <- distinct(input_df_RM)
  df_RM$kRM <- exp(predict(mod_k_RM, newdata=df_RM, re.form = ~0))
  
  output_kM <- left_join(df_AM, df_RM, by = c("date", "tempC", "Qout_ls"))%>%
    dplyr::select(tempC, Qout_ls, kRM, kAM)
  
  return(output_kM)
}