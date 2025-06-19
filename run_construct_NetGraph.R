# 0_load_packages ---------------------------------------------------------
library(brainGraph) # compute network properties
# remotes::install_version("igraph", version = "1.6.0")
library(igraph) # 1.6.0
packageVersion("igraph") # 1.6.0
library(data.table) # data.table
library(parallel)
library(doMC)
registerDoMC(14)

options(bg.subject_id='participant_id', bg.group='all_group')
grps = c('health_before', 'patient_before', 'health_after', 'patient_after')
# 1_network_construction --------------------------------------------------

## Network contains positive & negative & absolute
densities <- seq(0.25, 0.30, 0.05)

"my_data = fread('brainnetome_surface.csv')
my_data$x.mni <- as.numeric(my_data$x.mni)
my_data$y.mni <- as.numeric(my_data$y.mni)
my_data$z.mni <- as.numeric(my_data$z.mni)
my_data$lobe <- as.factor(my_data$lobe)
my_data$hemi <- as.factor(my_data$hemi)
my_data$gyrus <- as.factor(my_data$gyrus)
my_data$Yeo_7network <- as.factor(my_data$Yeo_7network)
my_data$Yeo_17network <- as.factor(my_data$Yeo_17network)
setkey(my_data, index)
setindex(my_data, name)
my_data
brainnetome_surface = my_data
save(brainnetome_surface, file='brainnetome_surface.rda')"

load('brainnetome_surface.rda')
brainnetome_surface
## --------------------------------------------------------------------------
save2dataframe <- function(data, file_path_suffix){
  vertex_df <- rbindlist(lapply(data, vertex_attr_dt))
  graph_df <- rbindlist(lapply(data, graph_attr_dt))
  fwrite(vertex_df, paste0(file_path_suffix, '_vertex.csv'), sep=',', col.names=TRUE)
  fwrite(graph_df, paste0(file_path_suffix, '_graph.csv'), sep=',', col.names=TRUE)
}
#### CFS: functional network graph weighted subject level
covars.all <- fread('all_subs_152_info.csv') # subs=152
inds = lapply(grps, function(x) covars.all[all_group == x, which = TRUE]) # 1=health=39, 2=patient=40
matfiles <- paste0('NetResults_invnodal/', covars.all$participant_id, '/raw_func_', covars.all$participant_id, '_invnodal.txt') 
my.mats <- create_mats(matfiles, modality = 'fmri',threshold.by = 'density',
                       mat.thresh = densities, inds = inds)
gw.sub <- vector('list', length(densities)) # ws: weighted subject
for (i in seq_along(densities)){
  gw.sub[[i]] <- make_brainGraphList(my.mats$A.norm.sub[[i]], 'brainnetome_surface' , level='subject',
                                     modality = 'fmri',threshold = densities[i],
                                     weighted = TRUE, gnames = covars.all$participant_id,
                                     grpNames = covars.all$all_group )
}
saveRDS(gw.sub, file=file.path('BrainGraphRDS/', 'CFS_func_weighted_subject.rds'), compress = 'xz')
save2dataframe(gw.sub, file=file.path('BrainGraphRDS/', 'CFS_func_weighted_subject'))

#### CFS: fiber network graph weighted subject level
matfiles <- paste0('NetResults_invnodal/', covars.all$participant_id, '/raw_fiber_', covars.all$participant_id, '_invnodal.txt') 
my.mats <- create_mats(matfiles, modality = 'dti',threshold.by = 'density',
                       mat.thresh = densities, inds = inds)
gw.sub <- vector('list', length(densities)) # ws: weighted subject
for (i in seq_along(densities)){
  gw.sub[[i]] <- make_brainGraphList(my.mats$A.norm.sub[[i]], 'brainnetome_surface', level='subject',
                                     modality = 'dti',threshold = densities[i],
                                     weighted = TRUE, gnames = covars.all$participant_id,
                                     grpNames = covars.all$all_group )
}
saveRDS(gw.sub, file=file.path('BrainGraphRDS/', 'CFS_fiber_weighted_subject.rds'), compress = 'xz')
save2dataframe(gw.sub, file=file.path('BrainGraphRDS/', 'CFS_fiber_weighted_subject'))

## --------------------------------------------------------------------------

########## end. author@kangwu