## Functions used in the Network-MB analysis script
## Created based on structure created by LE Koenig 2019
## Danielle Hare UConn
## Last Updated June 2022

## add in details here 

#########################
## Create Network Direction 
#########################

netset <- function(network, bq_m3hrm){
  # #Watershed Parameters Wolheim 2006
  # a=7.3 #NC, Heton 2011 
  # b=0.45 #NC, Helton 2011
  # c=0.408 #Wollheim 
  # d=0.294
  
  # Define hydrologic inflows/outflows for each reach (m3 hr-1 per m) 
  V(network)$m <- makeVertexAtt(network, bq_m3hrm, "estimate", by.df = "basin_id", by.g = "basin_id")
  V(network)$Qlocal <-  V(network)$length_reach * V(network)$m
  V(network)$Qest <- V(network)$lengthup_m * V(network)$m
  
# move through network order matters
for(i in 1:length(V(network))){
    #all values contributing to point
    
  up <- tryCatch(
    {as.numeric(unlist(strsplit(V(network)[i]$up, split= ",")))
      
    },
    error=function(cond) {
      return(NA)
    })
    
    # Discharge inflow from upstream network (m3 hr-1):
    if(is.na(up) == FALSE){
      V(network)$Qup[i] <- sum(V(network)$Qout[up]) #only one or junction will have two
    }  else{
      V(network)$Qup[i] =  0 #for origin nodes
    }
      
    
    #Discharge outflow to downstream reach (m3 hr-1):

    V(network)$Qout[i] <- sum(V(network)$Qlocal[i], V(network)$Qup[i], na.rm = TRUE)#if_else(V(network)$Qlocal[i]> 0, sum(V(network)$Qlocal[i], V(network)$Qup[i], na.rm = T), V(network)$Qup[i]) #as long as there is a postive Qlocal, add Qlocal + Qup to get Qout, otherwise just use Qup - accounts for losing reaches 
    # V(network)$width_m[i] <- if_else(V(network)$Qout[i]> 0, a * (V(network)$Qout[i]/3600)^b, 0) #already in m3/hr, needs to be in m3/s
    # V(network)$depth_m[i] <- if_else(V(network)$Qout[i]> 0, c * (V(network)$Qout[i]/3600)^d, 0)
    V(network)$Bedarea_m2[i] <- V(network)$width_m[i] * V(network)$length_reach[i]
    V(network)$CSarea_m2[i] <- V(network)$width_m[i] * V(network)$depth_m[i]
    V(network)$avg_v[i] <- if_else(V(network)$Qout[i]> 0, V(network)$Qout[i] / V(network)$CSarea_m2[i], 0)
    #V(network)$runoff_mday[i] = V(network)$Qout[i]/V(network)$CatchmentA[i] #need to check q_ma versus q0001e
 }
 
  return(network)
}

#add up stream node values to a network and calculate reach length into a node based on edges in
up_nodes <- function(network){
  up_list <- lapply(1:length(V(network)),function(i){
    #all values contributing to point
    up <- igraph::neighbors(network, i,  mode=c("in"))

    return(up)
    
  })
  
  in_edge <- lapply(1:length(V(network)),function(i){
    e_ID <- incident(network, i, mode = c("in"))
    return(sum(e_ID$length_m))
  })
  #
  V(network)$up <- sapply(up_list, function(s) if (length(s) == 0) NA_character_ else paste(s, collapse = ","))
  V(network)$length_reach <- sapply(in_edge, function(s) if (s == 0) 10 else as.numeric(paste(s, collapse = ","))) #if there is no edges in make the length in 10m
  
  return(network)
  
}


## ---------- Functions used in the mass balance model --------- ##
move_OC <- function(network, network_pre){
  V(network)$ind_order <- seq(1, length(V(network)))
  V(network)$FTOC_up_A <- 0
  V(network)$FTOC_up_A <- 0
  V(network)$DOC_up <- 0
  V(network)$FTOC_out_R <- 0
  V(network)$FTOC_out_R <- 0
  V(network)$DOC_out <- 0
  V(network)$FTOC_out <- 0
  
  
  #move through nodes - this network is structured so nodes move from headwaters to outlet, thus has to be a loop (order matters)
for(i in seq_along(V(network))){
  #for(i in 1:length(V(network))){
#print(V(network)$ind_order[i])
#print(V(network)$name[i])
    ##NOTE FPOC is in ADFM and DOC is in gm3
    if(V(network)[i]$up.all > 1){ #if there are contributing nodes (not orgin sites)
      up <- as.numeric(unlist(strsplit(V(network)[i]$up, split= ","))) #directly contributing nodes - assumes every timestep (here hourly, all transported carbon moves to next downstream reach)
      V(network)$DOC_up[i] <- sum(V(network_pre)$DOC_out[up], na.rm = TRUE) #sum the DOC_out from the previous timestep - amount that has come to the reach
      V(network)$DOC_out[i] <- (V(network)$DOC_up[i] + V(network)$DOC_local_gC[i]) # * (1-exp(-V(network)$vf[i]/V(network)$HL[i])) #add DOC from upstream from previous timestep to the DOC local from the current timestep
      # Exported carbon load to downstream reach (g d-1):
      V(network)$FTOC_up_A[i] <- sum(V(network_pre)$FTOC_out_A[up], na.rm = TRUE)
      V(network)$FTOC_out_A[i] <- V(network)$FTOC_up_A[i] + V(network)$FTOC_local_A[i]
      V(network)$FTOC_up_R[i] <- sum(V(network_pre)$FTOC_out_R[up], na.rm = TRUE)
      V(network)$FTOC_out_R[i] <- V(network)$FTOC_up_R[i] + V(network)$FTOC_local_R[i]
      V(network)$FTOC_out[i] <- V(network)$FTOC_out_R[i] + V(network)$FTOC_out_A[i]
    }  else{#origin sites - no contributing streams
      V(network)$DOC_out[i] <- V(network)$DOC_local_gC[i]
      V(network)$FTOC_out_A[i] <- V(network)$FTOC_local_A[i]
      V(network)$FTOC_out_R[i] <- V(network)$FTOC_local_R[i]
      V(network)$FTOC_out[i] <- V(network)$FTOC_local_A[i] + V(network)$FTOC_local_R[i]
    }
  }

  return(network)
}


### previous time step
previous_ts <- function(network) {
      p_network <- NULL
      function() {
        p_network <<- network
      }
      p_network
    }

## ----------- Calcuate radian date from julian day------ ##
rad_day <- function(x){
  rad_day1 <- 2*pi*x/365
  return(rad_day1)
}


#--------- Summarize the ts_all file to look at Watershed Carbon Fluxes
daily_clean_summarise <- function(ts_all){
  df <-  ts_all %>%
    group_by(name, date) %>% #summarize by node and date
    summarise(
      tempC = mean(tempC),
      Q_outlet = max(Qout),
      kA_avg = mean(kAt),
      kR_avg = mean(kRt),
      C_StStock_gC_avg = mean(POC_sStock_AFDMg_slow + POC_sStock_AFDMg_fast, na.rm = TRUE) * 0.484, #take mean Standing Stock for each node for each day
      C_LitterIn_gC = sum(ClocalLit_AFDMg_slow + ClocalLit_AFDMg_fast, na.rm = TRUE) * 0.484,
      C_breakdownTQ_gC = sum(POC_loss_g_C_At + POC_loss_gC_Rt, na.rm = TRUE),
      C_breakdownTQ_gC_frag = sum(POC_loss_gC_AF + POC_loss_gC_RF, na.rm = TRUE),
      C_breakdownTQ_gC_micrb = sum(POC_loss_gC_AM + POC_loss_gC_RM, na.rm = TRUE),
      C_gw_gC = sum(DOC_local_gC, na.rm = TRUE),
      WS_bentic_area_m2_avg =  mean(Bedarea_m2, na.rm = TRUE)
    )%>%
    ungroup() %>%
    group_by(date) %>%
    summarise(
      tempC = mean(tempC),
      Q_outlet = max(Q_outlet),
      kA_avg = mean(kA_avg),
      kR_avg = mean(kR_avg),
      C_StStock_gC_total = sum(C_StStock_gC_avg), #Sum the avg daily standing stock for the whole watershed
      C_LitterIn_gC = sum(C_LitterIn_gC),
      C_breakdownTQ_gC = sum(C_breakdownTQ_gC, na.rm = TRUE),
      C_breakdownTQ_gC_frag = sum(C_breakdownTQ_gC_frag, na.rm = TRUE),
      C_breakdownTQ_gC_micrb = sum(C_breakdownTQ_gC_micrb, na.rm = TRUE),
      C_gw_gC = sum(C_gw_gC, na.rm = TRUE),
      WS_bentic_area_m2 =  sum(WS_bentic_area_m2_avg, na.rm = TRUE)
    )
  
}

##keep for reference
# ##netset_replaced <- function(network, bq_m3dm){
#   #Watershed Parameters Wolheim 2006
#   a=7.3 #NC, Heton 2011 
#   b=0.45 #NC, Helton 2011
#   c=0.408 
#   d=0.294
#   
#   for(i in 1:length(V(network))){
#     
#     # Calculate mass-balance for each reach moving down the network from headwaters to mouth:
#     n_ID <- V(network)$name[i] #node ID being run
#     e_ID <- incident(network, n_ID, mode = c("in")) #get edges so the data is shared - this should replace the mutate to nodes. 
#     
#     #all values contributing to point
#     up <- igraph::neighbors(network, i,  mode=c("in")) #only single up (will be mostly 2)
#     ##put as a vertice attribute, want a string within a single column
#     #V(network)$upV[i] <- c(up)
#     #up.all.nodes <- head(unlist(ego(network,order=length(V(network)),nodes=V(network)[i],mode=c("in"),mindist=0)),-1)
#     
#     #e_ID$SEGMENT_ID
#     length_reach <- if(length(e_ID)>0){
#       sum(e_ID$length_m)} else(length_reach = 10)
#     #ifelse(length_reach == 0, 10, length_reach) # locations within any edges in are springs and thus have 10 meter contributing
#     V(network)$length_reach[i] <- length_reach
#     
#     # Define hydrologic inflows/outflows for each reach (m3 d-1)
#     # Discharge inflow from local catchment (m3 d-1):
#     #V(network)$Qlocal[i] <- V(network)$runoff_mday[i] * (V(network)$areasqkm[i]*10^6)
#     # #filter based on basin
#     V(network)$basin_id[i] <- if_else(V(network)$basin_id[i] == "CoweetaCreek", "ShopeFork", V(network)$basin_id[i]) #right now putting shope projection on the base as that incompasses most of the base
#     # 
#     bq <- bq_m3hrm %>%
#       dplyr::filter(basin_id == V(network)$basin_id[i]) %>% #average Baseflow m3/hr
#       mutate(m = estimate)
#     #b = b
#     
#     V(network)$Qlocal[i] <-  length_reach * bq$m #+ bq$b 
#     
#     
#     # Discharge inflow from upstream network (m3 hr-1):
#     if(length(up)>0){
#       V(network)$Qup[i] <- sum(V(network)$Qout[up]) #only one or junction will have two
#     }  
#     
#     #Discharge outflow to downstream reach (m3 hr-1):
#     V(network)$Qout[i] <- if_else(V(network)$Qlocal[i]> 0, sum(V(network)$Qlocal[i], V(network)$Qup[i], na.rm = T), V(network)$Qup[i])
#     V(network)$width_m[i] <- if_else(V(network)$Qout[i]> 0, a * (V(network)$Qout[i]/3600)^b, 0) #already in m3/hr, needs to be in m3/s
#     V(network)$depth_m[i] <- if_else(V(network)$Qout[i]> 0, c * (V(network)$Qout[i]/3600)^d, 0)
#     V(network)$Bedarea_m2[i] <- V(network)$width_m[i] * length_reach
#     V(network)$CSarea_m2[i] <- V(network)$width_m[i] * V(network)$depth_m[i]
#     V(network)$avg_v[i] <- if_else(V(network)$Qout[i]> 0, V(network)$Qout[i] / V(network)$CSarea_m2[i], 0)
#     #V(network)$runoff_mday[i] = V(network)$Qout[i]/V(network)$CatchmentA[i] #need to check q_ma versus q0001e
#   }
#   return(network)
# }