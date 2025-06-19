# 0_load_packages ---------------------------------------------------------
library(brainGraph) # compute network properties
# remotes::install_version("igraph", version = "1.6.0")
library(igraph) # 1.6.0
packageVersion("igraph") # 1.6.0
library(data.table) # data.table
library(parallel)
library(doMC)
registerDoMC(10)

options(bg.subject_id='participant_id', bg.group='group')
grps = c('health', 'patient')
# 1_network_construction --------------------------------------------------

## Network contains positive & negative & absolute
densities <- seq(0.10, 0.35, 0.05)
atlas = 'brainnetome'

matfiles <- Sys.glob('NetResults_invnodal/*/raw_func_*_invnodal.txt')
thr_25_list = vector()
for (i in matfiles){
  All_value = fread(i, header = F) %>%
    as.matrix(.) %>%
    get_thresholds(., densities)
  thr_25 = All_value[4]
  output = paste(i, thr_25)
  print(output)
  thr_25_list[i] = thr_25
}

max(thr_25_list)
min(thr_25_list)
