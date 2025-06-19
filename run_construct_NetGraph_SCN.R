library(data.table)
library(brainGraph)

OS <- .Platform$OS.type
if (OS == 'windows'){
  pacman::p_load(snow, doSNOW)
  num.cores <- as.numeric(Sys.getenv('NUMBER_OF_PROCESSORS'))
  cl <- makeCluster(num.cores, type='SOCK')
  clusterExport(cl, 'sim.rand.graph.par')
  registerDoSNOW(cl)
} else {
  suppressMessages(library(doMC))
  registerDoMC(detectCores()-1)
}

options(bg.subject_id='participant_id', bg.group='all_group')
lhrh <- fread('BrainGraphRDS/CFS_thickness_mat_exclude_Subject_019.csv')
# sort name
load('brainnetome_surface.rda')
assign("brainnetome", brainnetome_surface, envir = .GlobalEnv)
setkey(brainnetome_surface, index)
newname <- c('participant_id', brainnetome_surface$name)
colnames(lhrh) <- newname
# begin
covars.all <- fread('subs_146_info.csv')
covars.all = covars.all[subject!='Subject_019', ]
covars.all[, sex := as.factor(sex)]
covars <- covars.all[,c('age', 'sex', 'BMI', 'participant_id', 'all_group')]

myResids <- get.resid(lhrh, covars, atlas = 'brainnetome', exclude.cov = 'all_group')

densities <- seq(0.01, 0.34, 0.01)
corrs <- corr.matrix(myResids, densities=densities)

g <- lapply(seq_along(densities), function(x) 
  make_brainGraphList(corrs[x], modality='thickness'))

dt.G <- rbindlist(lapply(g, graph_attr_dt))  
dt.V <- rbindlist(lapply(g, vertex_attr_dt))
fwrite(dt.G, 'BrainGraphRDS/CFS_thickness_scn_graph.csv')
fwrite(dt.V, 'BrainGraphRDS/CFS_thickness_scn_vertex.csv')

########################################
options(bg.subject_id='participant_id', bg.group='all_group')
lhrh <- fread('BrainGraphRDS/CFS_coupling_mat_exclude_Subject_019.csv')
# sort name
load('brainnetome_surface.rda')
assign("brainnetome", brainnetome_surface, envir = .GlobalEnv)
setkey(brainnetome_surface, index)
newname <- c('participant_id', brainnetome_surface$name)
colnames(lhrh) <- newname

# begin
covars.all <- fread('subs_146_info.csv')
covars.all = covars.all[subject!='Subject_019', ]
covars.all[, sex := as.factor(sex)]
covars <- covars.all[,c('age', 'sex', 'BMI', 'participant_id', 'all_group')]

myResids <- get.resid(lhrh, covars, atlas = 'brainnetome', exclude.cov = 'all_group')

densities <- seq(0.01, 0.34, 0.01)
corrs <- corr.matrix(myResids, densities=densities, what = 'raw')  # raw

g <- lapply(seq_along(densities), function(x) 
  make_brainGraphList(corrs[x], modality='coupling'))

dt.G <- rbindlist(lapply(g, graph_attr_dt))  
dt.V <- rbindlist(lapply(g, vertex_attr_dt))
fwrite(dt.G, 'BrainGraphRDS/CFS_coupling_scn_graph.csv')
fwrite(dt.V, 'BrainGraphRDS/CFS_coupling_scn_vertex.csv')

########################################
options(bg.subject_id='participant_id', bg.group='all_group')
lhrh <- fread('BrainGraphRDS/CFS_coupling_mat_exclude_Subject_019.csv')
# sort name
load('brainnetome_surface.rda')
assign("brainnetome", brainnetome_surface, envir = .GlobalEnv)
setkey(brainnetome_surface, index)
newname <- c('participant_id', brainnetome_surface$name)
colnames(lhrh) <- newname

# begin
covars.all <- fread('subs_146_info.csv')
covars.all = covars.all[subject!='Subject_019', ]
covars.all[, sex := as.factor(sex)]
covars <- covars.all[,c('age', 'sex', 'BMI', 'participant_id', 'all_group')]

myResids <- get.resid(lhrh, covars, atlas = 'brainnetome', exclude.cov = 'all_group')

densities <- seq(0.01, 0.34, 0.01)
corrs <- corr.matrix(myResids, densities=densities)  # raw

g <- lapply(seq_along(densities), function(x) 
  make_brainGraphList(corrs[x], modality='coupling'))

dt.G <- rbindlist(lapply(g, graph_attr_dt))  
dt.V <- rbindlist(lapply(g, vertex_attr_dt))
fwrite(dt.G, 'BrainGraphRDS/CFS_coupling_residual_scn_graph.csv')
fwrite(dt.V, 'BrainGraphRDS/CFS_coupling_residual_scn_vertex.csv')

########## end. author@kangwu