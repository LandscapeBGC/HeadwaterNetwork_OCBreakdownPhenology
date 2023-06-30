
library(igraph)        # network analysis; create river network graphs
library(tidygraph)
library(tidyverse)
library(tidygraph)
library(broom)
library(ggnetwork)    #clean up igraph presentation?
library(ggplot2)
theme_set(theme_bw())
library(lubridate)
library(dplyr) 
library(lme4)
library(merTools)
# general data manipulation and cleaning
#library(magick) #create animation
#source("R/igraphNetwork_MB.R")
source("R/FUNCTIONS_MB.R")
source("R/FUNCTIONS_kPrediction.R")
library('devtools') 

library('Rsenal') ##install_github('brooksandrew/Rsenal') 
library(fuzzyjoin)

Boltz <- (8.617E-5)
