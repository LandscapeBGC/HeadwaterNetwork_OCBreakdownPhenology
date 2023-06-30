###### 
#Find avg flow for each section 
library(igraph)
library(tidyverse)
library(lubridate)
library(tidygraph)

net_lstQ <- readRDS("output/data/net_lst_baseQ_v2.RDS")

#make igraph vertices into one dataframe
net_Q_df <- do.call(rbind, lapply(net_lstQ, igraph::as_data_frame, what = "vertices")) %>%
  mutate(month = month(date))

#This is done in two steps to account for the non-calendar year start dates and get better annual avg Q
avg_Q <- net_Q_df %>%
  group_by(name, month)%>%
  summarise(avgQ = mean(Qout))%>%
  group_by(name)%>%
  summarise(avgQ = mean(avgQ))%>%
  mutate(
            avg_width = 7.3 * (avgQ/3600)^0.45, #from Helton 2011 and Leopold 1953
            avg_depth = 0.25 * (avgQ/3600)^0.25 #from Catalan 2022 (Eq2) and Leopold 1953
    )

#Scaling Metrics
#Watershed Parameters Wolheim 2006
# a=7.3 #NC, Heton 2011 
# b=0.45 #NC, Helton 2011
# c=0.408 #Wollheim 
# d=0.294

#add annual avg Q back into igraph networks 
net_Q_lst <- lapply(net_lstQ, function(x){
    tidygraph::as_tbl_graph(x) %>%
    activate(nodes) %>%
      dplyr::select(-"avgQ", -"avg_width", -"avg_depth")%>% #remove previous runs to replace with newest AvgQ values
    left_join(., y = avg_Q, by = "name")%>%
    mutate(
      width_m   = if_else(Qout> 0, avg_width * (((Qout/3600)/(avgQ/3600)))^0.26, 0), #from Catalan 2022 eq3
      depth_m   = if_else(Qout> 0, avg_depth * (((Qout/3600)/(avgQ/3600)))^0.4, 0), #from Catalan 2022 eq3
      Qout_ls   = Qout/3.6, #m3/hr to l/s
      # V(network)$width_m[i] <- if_else(V(network)$Qout[i]> 0, a * (V(network)$Qout[i]/3600)^b, 0) #already in m3/hr, needs to be in m3/s
      # V(network)$depth_m[i] <- if_else(V(network)$Qout[i]> 0, c * (V(network)$Qout[i]/3600)^d, 0)
      Bedarea_m2 = width_m * length_reach,
      CSarea_m2 = width_m * depth_m,
      avg_v_perh = if_else(Qout> 0, Qout/CSarea_m2, 0)
    )
}
)

net_lstQ <- write_rds(net_Q_lst, "output/data/net_lst_baseQ_v2.RDS")

review_df <- do.call(rbind, lapply(net_Q_lst, igraph::as_data_frame, what = "vertices")) %>%
  dplyr::select("stream", "n_lscp_name", "date", "Qout_ls")
saveRDS(review_df, "data/final_lmer_model_network/Qout_inputdata.RDS")


ws_36_width <- review_df %>%
  dplyr::filter(stream == "WS36")



