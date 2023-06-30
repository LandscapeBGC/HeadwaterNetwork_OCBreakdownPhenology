###################################
###################################
##Coweeta Carbon Network Analysis##
###################################
###################################
## Coweeta Average Air Temperature  12.8 C (cs01 met Jan 2003 - Dec 2018)

#source all req packages and scripts
source("R/global.R")
gc() #make sure all memory possible is available

#USER INPUT - Choose Litter Type Scenario
litter_type <- "Acer" #"Mixed", "Rhodo"

# Data created in NetworkOC_Supplemental_1.Rmd
if(litter_type == "Acer"){
  ####Direct and Lateral POC Inputs####
  #ClocalLit_AFDMg <- readRDS("data/ClocalLit_g.RDS")
  ClocalLit_AFDMg <- readRDS("data/output/POM_In.RDS") #%>%
  
  
  ##qual_In changed to 100% for 
  qual_IN <- readRDS("data/output/POM_In_speciesperc.RDS") %>%
    dplyr::select(1:2) %>%
    mutate(percent_low.fit = 0)#0 for acer, 1 for rhodo, comment out for mixed
  
  ClocalLit_AFDMg <- left_join(ClocalLit_AFDMg, qual_IN)%>%
    mutate(Cdirect_gm2hr_slow = Cdirect_gm2hr_ * percent_low.fit, #based on figure 2 Webster 2001
           Cdirect_gm2hr_fast = Cdirect_gm2hr_ * (1-percent_low.fit),
           Clateral_gmhr_slow = Clateral_gmhr_ * 0, #0.5 based on annual blow in from Webster 2001 *1 or *0 for acer alone or rhodo alone 05 for mixed
           Clateral_gmhr_fast = Clateral_gmhr_ * 1 #0.5 based on annual blow in from Webster 2001 *1 or *0 for acer alone or rhodo alone 05 for mixed
    )
  
}else({
  if(litter_type == "Rhodo"){
    ####Direct and Lateral POC Inputs####
    #ClocalLit_AFDMg <- readRDS("data/ClocalLit_g.RDS")
    ClocalLit_AFDMg <- readRDS("data/output/POM_In.RDS") #%>%
    
    
    ##qual_In changed to 100% for 
    qual_IN <- readRDS("data/output/POM_In_speciesperc.RDS") %>%
      dplyr::select(1:2) %>%
      mutate(percent_low.fit = 1)#0 for acer, 1 for rhodo, comment out for mixed
    
    ClocalLit_AFDMg <- left_join(ClocalLit_AFDMg, qual_IN)%>%
      mutate(Cdirect_gm2hr_slow = Cdirect_gm2hr_ * percent_low.fit, #based on figure 2 Webster 2001
             Cdirect_gm2hr_fast = Cdirect_gm2hr_ * (1-percent_low.fit),
             Clateral_gmhr_slow = Clateral_gmhr_ * 1, #0.5 based on annual blow in from Webster 2001 *1 or *0 for acer alone or rhodo alone 05 for mixed
             Clateral_gmhr_fast = Clateral_gmhr_ * 0 #0.5 based on annual blow in from Webster 2001 *1 or *0 for acer alone or rhodo alone 05 for mixed
      )
    
  }else({
    ####Direct and Lateral POC Inputs####
    ClocalLit_AFDMg <- readRDS("data/output/POM_In.RDS") #%>%

    qual_IN <- readRDS("data/output/POM_In_speciesperc.RDS") %>%
      dplyr::select(1:2)
    
    ClocalLit_AFDMg <- left_join(ClocalLit_AFDMg, qual_IN)%>%
      mutate(Cdirect_gm2hr_slow = Cdirect_gm2hr_ * percent_low.fit, #based on figure 2 Webster 2001
             Cdirect_gm2hr_fast = Cdirect_gm2hr_ * (1-percent_low.fit),
             Clateral_gmhr_slow = Clateral_gmhr_ * 0.5, #based on annual blow in from Webster 2001 *1 or *0 for acer alone or rhodo alone 05 for mixed
             Clateral_gmhr_fast = Clateral_gmhr_ * 0.5 #0.5 based on annual blow in from Webster 2001 *1 or *0 for acer alone or rhodo alone 05 for mixed
      )}
  )}
)


#Set start date
intial_dates = as.POSIXct(c("08-01-2018"), format = "%m-%d-%Y", tz = "GMT")

###########################################################################
## Run Multiple Temperature Scenarios in sequence

############################
#define hour time steps and units
##############################
timesteps = 2 * (24 * 365) -1 #minus one hour to end on correct day
ts_units = "hour" #hour'
################
#Breakdown Model
################
# Created within NetworkOC_Supplemental_4
# Scale values for TQ data for respective lmer model
scaled_TQ_df <- readRDS("data/output/scaled_TQ_df.RDS")
mod_k_AM <- readRDS("data/output/acer_microbes_model.RDS")
mod_k_AF <- readRDS("data/output/acer_shredder_model.RDS")
mod_k_RM <- readRDS("data/output/rhodo_microbes_model.RDS")
mod_k_RF <- readRDS("data/output/rhodo_shredder_model.RDS")

### Set up Temperature Scnearios
#Climate scenarios derived from Hare et al. 2021; 0.04 C/year (which is the mean rate of increase for both shall and atm), verus for 0.01 C/year; therefore for 50 years 2 C for atmospheric and shallow, versus 0.5 for deep groundwater 
# Created within NetworkOC_Supplemental_3
temp_sin <- readRDS("data/output/temp_sin.RDS")
scen_T <- list(
               base = temp_sin, #Observed Coweeta Stream Temperature
               deepGW = mutate(temp_sin, amp = 4, phase = 200, ymean = 12), #deep GW
               shalGW = mutate(temp_sin, amp = 5.5, phase = 220, ymean = 12),#shallow GW
               low_GW = mutate(temp_sin, amp = 9, phase = 200, ymean = 13),#minimal GW influence
               low_GW_2 = mutate(temp_sin, amp = 9, phase = 200, ymean = 15),#, air-coupled warming scenario
               deepGW_05 = mutate(temp_sin, amp = 4, phase = 200, ymean = 12.5), #deep GW warming scenario
               shalGW_2 = mutate(temp_sin, amp = 5.5, phase = 220, ymean = 14)#shallow GW warming scenario
               )
scen_temp <- names(scen_T) #list the scen
intial_objects <- ls() #create list of objects needed preRun, to clean up everything but whatis listed above on iterative scenario runs 

###############################
## Run each Thermal scenario defined in scen_T
###############################
for(scen in scen_temp){
  rm(list=setdiff(ls(), c(intial_objects, "scen", "intial_objects"))) #remove everything in runs before 
  
  
    modelrun_start <- now()# start timer
    #set up output dataframe
    network_ts_all <- list()
    network_ts_day <- list()
    network_ts_day_id <- list()
    
    ###################
    ##Input Variables## 
    ###################
    
    ##########################
    #### Read Input Data #####
    ##########################
    #Read in Daily Flow 
    net_lstQ <- readRDS("data/net_lst_baseQ_v2.RDS")#updated to reflect appropriate scaling, and is the network for each day
    temp_sin <- scen_T[scen][[1]]
    
    #Add daily Stream Temperature to Each Frame 
    net_lstQ <- lapply(net_lstQ, function(x){
      tidygraph::as_tbl_graph(x) %>%
        activate(nodes) %>%
        left_join(., temp_sin, by = c("n_lscp_name" = "stream"))%>%
        mutate(
          date = as.Date(date),
          day = yday(as.Date(date)),
          tempC = round(amp * cos(rad_day(day -  phase)) +  ymean,1),
          Qout_ls = round(Qout_ls, 3)
        )
    }
    )
    
    
    # ###read decomposition models
    ### Breakdown Inputs
    #Create Dataframe for Temp, Flow and Breakdown for two year model run (if model runs extends into 2024 has to be updated to include additional leap year)
    input_df_TQ <- do.call(rbind, lapply(net_lstQ , igraph::as_data_frame, what = "vertices")) %>% #calculating extactly tempC combos are within the network scenario run
      dplyr::select("stream", "n_lscp_name", "date", "Qout_ls", "tempC")%>%
      distinct()
      
    #calculate k model for temperature run
    message("Start Breakdown Calculation")
    k_frag <- fragmentation_k_calc(input_df_TQ, scaled_TQ_df, mod_k_AF, mod_k_RF)%>%
      dplyr::select("tempC", "Qout_ls", "kAF", "kAF_UCI", "kAF_LCI", "kRF", "kRF_UCI", "kRF_LCI")%>%
      distinct()
    message("Frag Breakdown Calculated")
    k_rates <- left_join(k_frag, microbial_k_calc(input_df_TQ, scaled_TQ_df, mod_k_AM, mod_k_RM))%>%
      distinct() #only calculate for unique temperature/discharge pairs - saves time
    rm(input_df_TQ, k_frag)
    message("End Breakdown Calculation")
    
    #CPOM # no longer used, but can create standing stock comparison 
    cpom_gm2 <- readRDS("data/cpom_gm2.RDS")
    #create a vector with the correct cpom initialization for each date. 
    intial_cpomSS = cpom_gm2 %>% dplyr::filter(Jdate %in% yday(intial_dates))
    
    #### Seep DOC Data by Month ####
    DOC_gw <- readRDS("data/output/DOC_seep_table.RDS")
#######################
        #set up run for fragmentation model mean, model UCI, and model LCI
frag_model <- c("UCI", "LCI", "M") 
        
        for(frag_model_i in frag_model){
    ####################################
    #SETTING UP INITIAL DATE FOR THE RUN,
              j = 1 #previously had been running multiple start dates, but now single run 
      
              s_date = as.POSIXct(intial_dates[j], format = "%m-%d-%Y %H:%M", tz = "GMT")#("10-01-2018 00:00", format = "%m-%d-%Y %H:%M") #define model start date
              s_Jday = yday(s_date) #starting Julian Day 
              s_month = month(s_date) #starting month
              
              ##for loop for each time step
              #water yield, doc yield, litter inputs, temperature
                      dates <- seq(from = as.POSIXct(s_date, tz = "GMT"), to = as.POSIXct(s_date + hours(timesteps), tz = "GMT"), by = "hour")
                      
                      #remove leap day
                      dates1 <- as.data.frame(dates)%>%
                        .[!(format(.$dates,"%m") == "02" & format(.$dates, "%d") == "29"), , drop = FALSE]
                      dates <- dates1$dates #make back into a vector
                      
                      #env_network <- new.env(parent = emptyenv())
                      
                      
              #REMOVE OBJECT TO REFRESH ENVIRONMENT        
                    if(exists("network_pre")){#, envir = env_network)){
                      
                      #rm("network_pre", envir = env_network)
                      rm(network_pre, inherits = TRUE) #If inherits is TRUE then parents of the supplied directory are searched until a variable with the given name is encountered.
                      
                      } #remove from env
                      

              ##HOURLY TIME STEPS        
              net_lst <- lapply(dates, function(t_s, env = parent.frame(), inherits = FALSE){#env = env_network, inherits = FALSE){#
                              #print(t_s)#print current time step
                              ##what was the previous time step , fifelse presevers type and class of inputs this catches daylight savings, but using GMT now (5-15-2023)
                              ts_pre <- dplyr::if_else(is.na(t_s - hours(1)) == TRUE, t_s - hours(2), t_s - hours(1)) 
                              
                              #initialize date details
                              #pull timestep specific values for filtering 
                              #year_run <- year(t_s)
                              #for longer model runs than 1 year, after july 31, 2019, change back to 2018 - this is only for the Q discharge data
                              int_date = as.POSIXct("07-31-2019 23:00", format = "%m-%d-%Y %H:%M", tz = 'GMT')#("10-01-2018 00:00", format = "%m-%d-%Y %H:%M")
                              leap_date = as.POSIXct(paste0("02-29-", year(t_s), " 23:00" , tz = 'GMT'), format = "%m-%d-%Y %H:%M")  
                              #t_s_Q <- if_else(t_s > int_date, `year<-`(t_s, 2018), t_s)
                              
                              #dealing with leap year, correct each leap year to read days as non-leap (feb 29 removed and dec 31 is day 364)
                              leap_year(t_s)
                              #date to pull Discharge Data From 
                              date_cQ <- ifelse(leap_year(t_s) == TRUE && t_s > leap_date, as.character(yday(as.Date(t_s))-1), as.character(yday(as.Date(t_s))))#make is march 1 otherwise date march 1 become feb 29!
                              day <- ifelse(leap_year(t_s) == TRUE && t_s > leap_date, yday(as.Date(t_s))-1, yday(as.Date(t_s)))#make is march 1 otherwise date march 1 become feb 29!
                              month <- month(t_s)
                              #season <- quarter(t_s, fiscal_start = 1) #season starting with jan
              
                              
          ### INITIAL INPUT
                    ### Test if there is a previous time step, or if its a new day. As this generates a new Q, temp, breakdown based on day
                    if(!exists("network_pre") || yday(t_s) != yday(ts_pre)){
                      
                      #print every month
                      if(day(as.Date(t_s)) == 1){
                      message(t_s)
                      message(yday(t_s))
                      message(now())
                      }
                      #######################
                      ### Discharge ######
                      ####################
                      #if there is no previous network (eg start of the session), or there is a new day pull the correct date baseQ (standing stock will be still use previous timestep)
                      network <- net_lstQ[[date_cQ]] #pull the correct network generated from baseQ data
                      
                      #add date to igraph
                      V(network)$date <- as.character(t_s)  #, format = "%Y-%m-%d")
                      V(network)$date_d <- as.character(as.Date(t_s))

                          #################
                          ### CPOM ###
                          ################
                              ###find day landscape cpom standingstock g per m2 - interpolated and averaged
                              cpom_ts <- cpom_gm2 %>%
                                dplyr::filter(Jdate == day)
                          
                            ### Add Litter POC (direct and lateral) input
                              POM_input_day <- ClocalLit_AFDMg %>%
                                dplyr::filter(Jdate == day)
                              
                             ### Add Terrestrial Input DOC
                              DOC_seep_table_mon <- DOC_gw %>%
                                filter(mon == month)
    
                                
                                #add CPOM standing stock based 
                                V(network)$ss_POC_l <- cpom_ts$cpom_fit   #makeVertexAtt(network, df=cpom_ts, vname='cbom')#vname='cbom.afdm.g.m2', by.df='stream', by.g='n_lscp_name')
                                V(network)$ss_POC <- V(network)$Bedarea_m2 * V(network)$ss_POC_l #initial POC standing stock, only used in first time step
                                
                                #POC in - using both direct and lateral is essential, as direct is only m2 which is strongly dependant on the width estimate, while lateral is consistent as assocaited wiht length
                                V(network)$ClocalLit_AFDMg_slow <- (POM_input_day$Cdirect_gm2hr_slow *V(network)$Bedarea_m2) + #direct
                                  (POM_input_day$Clateral_gmhr_slow * V(network)$length_reach * 2) 
                                
                                
                                V(network)$ClocalLit_AFDMg_fast <- (POM_input_day$Cdirect_gm2hr_fast *V(network)$Bedarea_m2) + #direct
                                  (POM_input_day$Clateral_gmhr_fast * V(network)$length_reach * 2)
                                
                                
                                #If the flow is not negative, mutliply DOC seep value by reach Q added to get DOC from GW for the reach
                                V(network)$DOC_local_gC <- if_else(DOC_seep_table_mon$doc_ppm_av * V(network)$Qlocal < 0, 0, DOC_seep_table_mon$doc_ppm_av * V(network)$Qlocal) # *1000 / 1000 mg / L <- g/m3
                                
                            
                                ###################
                                ## Breakdown #####
                                #################
                                ##k using carolyn's models##
                                ##updated 2023-04-26## 
                                
                                # Retrive Calculated Breakdown Rate 
                                # uses tidygraph
                                
                                network <- as_tbl_graph(network) %>%
                                  activate(nodes) %>%
                                  mutate(tempC = round(tempC, 1),
                                         Qout_ls = round(Qout_ls, 3))%>%
                                  left_join(., k_rates, by = c("tempC", "Qout_ls")) #dropped "'date_d' = "date" 
                                
                                        #################################
                              
                                #output into the network 
                                if(frag_model_i == "LCI"){
                                  V(network)$kRt <- V(network)$kRM + V(network)$kRF_LCI
                                }else
                                  if(frag_model_i == "UCI"){
                                    V(network)$kRt <- V(network)$kRM + V(network)$kRF_UCI
                                  }else{
                                    V(network)$kRt <- V(network)$kRM + V(network)$kRF
                                  }
                                
                                
                                #
                                if(frag_model_i == "LCI"){
                                  V(network)$kAt <- V(network)$kAM + V(network)$kAF_LCI
                                }else
                                  if(frag_model_i == "UCI"){
                                    V(network)$kAt <- V(network)$kAM + V(network)$kAF_UCI
                                  }else{
                                    V(network)$kAt <- V(network)$kAM + V(network)$kAF
                                  }
                                
                                
                                
                                #############################
                                ##Indicates this part of if statment run if there is no pre network, or the network is from the same day - no new Q, T or k values
                                #message("No PRE or New day") #for bebugging
                                
                                
                        } else{ #USE THE PREVEIOUS DATA- except STandingstock 
                                network <- network_pre
                                #message("Same day") #for debugging
                                
                        }
                              
                          #########################################
                          ##-----POC STANDING STOCK CALC-------- ##
                          #########################################
                              #POC IN hourly, sum direct input and lateral - assume direct for full width 
                              ### Previous Time Step to determine present standing stock for reach length (sStock)
                              if(!exists("network_pre")){
    
                                #have no breakdown first timestep so the intial value should be the same for all scenarios. 
                                #V(network)$POC_sStock_AFDMg <-  V(network)$ClocalLit_AFDMg  #V(network)$ss_POC #+ #no longer start with value as I run for 100 days. 
                                V(network)$POC_sStock_AFDMg_slow <-  V(network)$ClocalLit_AFDMg_slow
                                V(network)$POC_sStock_AFDMg_fast <-  V(network)$ClocalLit_AFDMg_fast
                                
                                #Acer "fast"
                                V(network)$POC_loss_AFDMg_AF <-  0
                                V(network)$POC_loss_gC_AF <- 0
                                V(network)$FTOC_local_A <- 0
                                V(network)$POC_loss_AFDMg_AM <-  0
                                V(network)$POC_loss_gC_AM <- 0
                                V(network)$POC_loss_AFDMg_At   <- 0
                                V(network)$POC_loss_g_C_At   <- 0
                                
                                #Rhodo "slow"
                                V(network)$POC_loss_AFDMg_RF <-  0
                                V(network)$POC_loss_gC_RF <- 0
                                V(network)$FTOC_local_R <- 0
                                V(network)$POC_loss_AFDMg_RM <-  0
                                V(network)$POC_loss_gC_RM <- 0
                                V(network)$POC_loss_AFDMg_Rt   <- 0
                                V(network)$POC_loss_gC_Rt   <- 0
                            
                              } else({
                                #use previous time step POC standingstock + litter in 
                                V(network)$POC_sStock_AFDMg_slow <-  V(network_pre)$POC_AFDMg_slow + V(network)$ClocalLit_AFDMg_slow
                                V(network)$POC_sStock_AFDMg_fast <-  V(network_pre)$POC_AFDMg_fast + V(network)$ClocalLit_AFDMg_fast
                                
                                ##################################################
                                #Acer "fast"
                                if(frag_model_i == "LCI"){
                                V(network)$POC_loss_AFDMg_AF <-  V(network)$POC_sStock_AFDMg_fast * (V(network)$kAF_LCI/24)
                                }else
                                  if(frag_model_i == "UCI"){
                                                            V(network)$POC_loss_AFDMg_AF <-  V(network)$POC_sStock_AFDMg_fast * (V(network)$kAF_UCI/24)
                                  } else{
                                          V(network)$POC_loss_AFDMg_AF <-  V(network)$POC_sStock_AFDMg_fast * (V(network)$kAF/24)
                                  }
                                V(network)$POC_loss_gC_AF <- V(network)$POC_loss_AFDMg_AF *0.484 #convert gC
                                V(network)$FTOC_local_A <- V(network)$POC_loss_AFDMg_AF # I have two for clarity right now when conisdering transport in move_OC
                                V(network)$POC_loss_AFDMg_AM <-  V(network)$POC_sStock_AFDMg_fast * (V(network)$kAM/24)
                                V(network)$POC_loss_gC_AM <- V(network)$POC_loss_AFDMg_AM *0.484 #convert to gC
                                V(network)$POC_loss_AFDMg_At   <- V(network)$POC_loss_AFDMg_AF + V(network)$POC_loss_AFDMg_AM
                                V(network)$POC_loss_g_C_At   <- V(network)$POC_loss_AFDMg_At *0.484#convert to gC
                                
                                ###########################
                                #Rhodo "slow"
                                if(frag_model_i == "LCI"){
                                  V(network)$POC_loss_AFDMg_RF <-  V(network)$POC_sStock_AFDMg_slow * (V(network)$kRF_LCI/24)
                                }else
                                  if(frag_model_i == "UCI"){
                                    V(network)$POC_loss_AFDMg_RF <-  V(network)$POC_sStock_AFDMg_slow * (V(network)$kRF_UCI/24)
                                  }else{
                                    V(network)$POC_loss_AFDMg_RF <-  V(network)$POC_sStock_AFDMg_slow * (V(network)$kRF/24)
                                  }
                                
                                V(network)$POC_loss_gC_RF <- V(network)$POC_loss_AFDMg_RF *0.484 #convert gC
                                V(network)$FTOC_local_R <- V(network)$POC_loss_AFDMg_RF # I have two for clarity right now when conisdering transport in move_OC
                                V(network)$POC_loss_AFDMg_RM <-  V(network)$POC_sStock_AFDMg_slow * (V(network)$kRM/24)
                                V(network)$POC_loss_gC_RM <- V(network)$POC_loss_AFDMg_RM *0.484 #convert to gC
                                V(network)$POC_loss_AFDMg_Rt   <- V(network)$POC_loss_AFDMg_RF + V(network)$POC_loss_AFDMg_RM
                                V(network)$POC_loss_gC_Rt   <- V(network)$POC_loss_AFDMg_Rt *0.484#convert to gC
                              
                              })
              
                              #POC standing Stock Loss per hour
                              V(network)$POC_AFDMg_slow <- V(network)$POC_sStock_AFDMg_slow - V(network)$POC_loss_AFDMg_Rt
                              V(network)$POC_AFDMg_fast <- V(network)$POC_sStock_AFDMg_fast - V(network)$POC_loss_AFDMg_At
                
                              ################
                              
              #   ### ADD IN FOR SERIAL###
              #   ####################################################
              #             ### Calculate movement within basin
              #             ### FPOC and DOC ###
              #             if(!exists("network_pre")){
              #               #set up transport
              #               V(network)$FTOC_up <- 0
              #               V(network)$DOC_up <- 0
              #               V(network)$ind_order <- seq(1, length(V(network)))
              #               V(network)$FTOC_up_A <- 0
              #               V(network)$FTOC_up_A <- 0
              #               V(network)$DOC_up <- 0
              #               V(network)$FTOC_out_R <- 0
              #               V(network)$FTOC_out <- V(network)$FTOC_local_A + V(network)$FTOC_local_R
              #               V(network)$DOC_out <- V(network)$DOC_local_gC
              # 
              #             } else {
              #                 network <- move_OC(network, network_pre)
              #                 #message("moveOC")
              #             }
              # ######################################################################
                              
                    ##set up environment for next timestep
                              
                              #YES ASSIGNING TO THE GLOBAL IS FROWNED ON _ WILL CHANGE TO DIFFERNT ONE BUT WORKS!!!!! 
                              assign("network_pre", network, envir = .GlobalEnv)#env_network
                              
                              if(year(t_s)>2018){
                                    network_df <- igraph::as_data_frame(network, what = "vertices")  
                                    # %>% 
                                    #                   dplyr::select(name, date_d, date, Qout, Bedarea_m2, 36:ncol(.))#make smaller for saving 
                                    # 
                                    # #message(Sys.time())
                                    return(network_df)
                              }else{
                                return(NULL)}
                              
                      })

              
              #remove network_pre for scenario running
            
              ##################
              ### OUTPUT #######
              ##################
              # #Create dataframe#
              
              #remoev first list in network list as there is no transport
              #net_lst[[1]] <- NULL
              
              #"convert igraphs to dataframe baed on vertices attribtues
              # ts_all <-  lapply(net_lst, function(i){ # start with second one as the first timestep doesnt have transport.
              #   igraph::as_data_frame(i, what = "vertices")
              # })
              message("Output Dataframes")
              net_lst <- net_lst[- c(seq(1:4320))]  #remove up until June 2019 for time and memory(2018 not saved at all)
              ts_all <-  do.call("rbind", net_lst)
                
              write_rds(ts_all, paste0("results/model_output/network_", litter_type, "_", scen, "_ts_all_",frag_model_i, "_nonserialC.RDS"))
              
              #summarize breakdown and gw_DOC
              ### Sum hourly timesteps to the day
              ts_day <- daily_clean_summarise(ts_all) 
              
              write_rds(ts_day, paste0("output_v2/data/network_", litter_type, "_", scen, "_ts_day_", frag_model_i, "_nonserialC.RDS"))
              ##Note the "lengthup_m" is from the original GIS calculations, so is determined seperately from this code. 
              ##Therefore, the SLout and lengthup_m should be ~ the same. 
              # #
              # show_env <- solveMB(){
              #   list(ran.in = environment(), 
              #        parent = parent.env(environment()), 
              #        objects = ls.str(environment()))
              # }
rm(ts_day, ts_all, network_pre, net_lst)
gc()
message(now() - modelrun_start)
        }# end for loop for CI levels 

}#end scenario lapply

