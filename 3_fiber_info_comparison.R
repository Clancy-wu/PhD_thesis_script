library(data.table)
options(warn = -1) 

## clinical info
clinic_info_73 = fread('subs_146_info.csv')
## fiber net
fiber_graph_file = 'BrainGraphRDS/CFS_fiber_weighted_subject_graph.csv'
fiber_graph = fread(fiber_graph_file)
fiber_graph = fiber_graph[threshold==0.25, ]
# [1] "threshold"            "Cp"                   "Lp"                   "E.global"             "mod"                 
# [6] "density"              "max.comp"             "num.tri"              "diameter"             "transitivity"        
# [11] "assort"               "strength"             "E.local.wt"           "E.global.wt"          "diameter.wt"         
# [16] "Lp.wt"                "num.hubs.wt"          "mod.wt"               "assort.lobe"          "assort.lobe.hemi"    
# [21] "asymm"                "spatial.dist"         "assort.gyrus"         "assort.Yeo_7network"  "assort.Yeo_17network"
# [26] "num.hubs"             "E.local"              "vulnerability"        "participant_id"       "atlas"               

select_col = 'E.global.wt'
net_df = fiber_graph
corr_between_net_behavior <- function(net_df, select_col){
  net_select_col = c(select_col, 'participant_id')
  net_select_df = net_df[, ..net_select_col]
  merge_df = merge(clinic_info_73, net_select_df, by='participant_id', all.x = T)
  setkey(merge_df, 'subject')
  merge_df_patient_before = merge_df[group=='patient'&treatment=='before', ]
  merge_df_patient_after = merge_df[group=='patient'&treatment=='after', ]
  merge_df_health_before = merge_df[group=='health'&treatment=='before', ]
  merge_df_health_after = merge_df[group=='health'&treatment=='after', ]
  ## compare
  p_pearson <- cor.test(
    merge_df_patient_after[,..select_col][[1]] - merge_df_patient_before[,..select_col][[1]],
    merge_df_patient_after[,.(FS14)][[1]] - merge_df_patient_before[,.(FS14)][[1]],
    method = 'pearson')
  p_spearman <- cor.test(
    merge_df_patient_after[,..select_col][[1]] - merge_df_patient_before[,..select_col][[1]],
    merge_df_patient_after[,.(FS14)][[1]] - merge_df_patient_before[,.(FS14)][[1]],
    method = 'spearman')
  
  h_pearson <- cor.test(
    merge_df_health_after[,..select_col][[1]] - merge_df_health_before[,..select_col][[1]],
    merge_df_health_after[,.(FS14)][[1]] - merge_df_health_before[,.(FS14)][[1]],
    method = 'pearson')
  h_spearman <- cor.test(
    merge_df_health_after[,..select_col][[1]] - merge_df_health_before[,..select_col][[1]],
    merge_df_health_after[,.(FS14)][[1]] - merge_df_health_before[,.(FS14)][[1]],
    method = 'spearman')
  
  p_h_t <- t.test(
    merge_df_patient_before[,..select_col][[1]], merge_df_health_before[,..select_col][[1]],
    paired = F, var.equal = T)
  p_h_w <- wilcox.test(
    merge_df_patient_before[,..select_col][[1]], merge_df_health_before[,..select_col][[1]],
    paired = F )  
  
  output <- paste('Changed Corr: \n    Patient: Pearsonr:', round(p_pearson$p.value,4), 
                  'Spearman:', round(p_spearman$p.value,4), 
                  '\n    Health: Pearsonr:', round(h_pearson$p.value,4), 
                  'Spearman:', round(h_spearman$p.value,4), 
                  '\nBefore Diff: \n    T test:', round(p_h_t$p.value,4),
                  'W test:', round(p_h_w$p.value,4), sep=' ' )
  return(cat(output))
}

corr_between_net_behavior(net_df, select_col)
## On fibertional network, no indicators cater with (1) changes corr with behavior, (2) baseline different

############################## vertex
clinic_info_73 = fread('subs_146_info.csv')
load('brainnetome_surface.rda')

fiber_vertex_file = 'BrainGraphRDS/CFS_fiber_weighted_subject_vertex.csv'
fiber_vertex = fread(fiber_vertex_file)
fiber_vertex = fiber_vertex[threshold==0.25, ]
# [1] "density"         "region"          "lobe"            "hemi"            "gyrus"           "Yeo_7network"   
# [7] "Yeo_17network"   "degree"          "comp"            "strength"        "knn.wt"          "s.core"         
# [13] "transitivity.wt" "E.local.wt"      "E.nodal.wt"      "Lp.wt"           "hubs.wt"         "comm.wt"        
# [19] "GC.wt"           "PC.wt"           "z.score.wt"      "asymm"           "dist"            "dist.strength"  
# [25] "knn"             "Lp"              "btwn.cent"       "hubs"            "ev.cent"         "lev.cent"       
# [31] "k.core"          "transitivity"    "E.local"         "E.nodal"         "vulnerability"   "eccentricity"   
# [37] "comm"            "GC"              "PC"              "z.score"         "participant_id"  "atlas"          
# [43] "modality"        "all_group"       "threshold"  

select_col = 'E.nodal'
net_df = fiber_vertex
corr_between_net_behavior_vertex <- function(net_df, select_col){
  net_select_col = c(select_col, 'region', 'participant_id')
  net_select_df = net_df[, ..net_select_col]
  merge_df = merge(clinic_info_73, net_select_df, by='participant_id', all.x = T)
  setkey(merge_df, 'subject')
  merge_df_patient_before = merge_df[group=='patient'&treatment=='before', ]
  merge_df_patient_after = merge_df[group=='patient'&treatment=='after', ]
  merge_df_health_before = merge_df[group=='health'&treatment=='before', ]
  merge_df_health_after = merge_df[group=='health'&treatment=='after', ]
  # all regions
  all_regions = brainnetome_surface$name
  region_p_pearson = vector()
  region_p_spearman = vector()
  region_h_pearson = vector()
  region_h_spearman = vector()
  region_p_h_t = vector()
  region_p_h_w = vector()
  
  for (i in seq_along(all_regions)){
    ## compare
    p_pearson <- cor.test(
      merge_df_patient_after[region==all_regions[i],..select_col][[1]] - merge_df_patient_before[region==all_regions[i],..select_col][[1]],
      merge_df_patient_after[region==all_regions[i],.(FS14)][[1]] - merge_df_patient_before[region==all_regions[i],.(FS14)][[1]],
      method = 'pearson')
    p_spearman <- cor.test(
      merge_df_patient_after[region==all_regions[i],..select_col][[1]] - merge_df_patient_before[region==all_regions[i],..select_col][[1]],
      merge_df_patient_after[region==all_regions[i],.(FS14)][[1]] - merge_df_patient_before[region==all_regions[i],.(FS14)][[1]],
      method = 'spearman')
    
    h_pearson <- cor.test(
      merge_df_health_after[region==all_regions[i],..select_col][[1]] - merge_df_health_before[region==all_regions[i],..select_col][[1]],
      merge_df_health_after[region==all_regions[i],.(FS14)][[1]] - merge_df_health_before[region==all_regions[i],.(FS14)][[1]],
      method = 'pearson')
    h_spearman <- cor.test(
      merge_df_health_after[region==all_regions[i],..select_col][[1]] - merge_df_health_before[region==all_regions[i],..select_col][[1]],
      merge_df_health_after[region==all_regions[i],.(FS14)][[1]] - merge_df_health_before[region==all_regions[i],.(FS14)][[1]],
      method = 'spearman')
    
    p_h_t <- t.test(
      merge_df_patient_before[region==all_regions[i],..select_col][[1]], merge_df_health_before[region==all_regions[i],..select_col][[1]],
      paired = F, var.equal = T)
    p_h_w <- wilcox.test(
      merge_df_patient_before[region==all_regions[i],..select_col][[1]], merge_df_health_before[region==all_regions[i],..select_col][[1]],
      paired = F ) 
    # save p value
    region_p_pearson[i] = p_pearson$p.value
    region_p_spearman[i] = p_spearman$p.value
    region_h_pearson[i] = h_pearson$p.value
    region_h_spearman[i] = h_spearman$p.value
    region_p_h_t[i] = p_h_t$p.value
    region_p_h_w[i] = p_h_w$p.value
  }
  # fdr
  #region_p_pearson_adj = p.adjust(region_p_pearson)
  #region_p_spearman_adj = p.adjust(region_p_spearman)
  #region_h_pearson_adj = p.adjust(region_h_pearson)
  #region_h_spearman_adj = p.adjust(region_h_spearman)
  #region_p_h_t_adj = p.adjust(region_p_h_t)
  #region_p_h_w_adj = p.adjust(region_p_h_w)
  out_df = data.table(
    region_p_pearson = region_p_pearson,
    region_p_spearman = region_p_spearman,
    region_h_pearson = region_h_pearson,
    region_h_spearman = region_h_spearman,
    region_p_h_t = region_p_h_t,
    region_p_h_w = region_p_h_w
  )
  return(out_df)
}

my_df = corr_between_net_behavior_vertex(net_df, select_col)
fwrite(my_df, file='E.nodal.csv')

##########################
clinic_info_73 = fread('subs_146_info.csv')
load('brainnetome_surface.rda')
Mean_Sd <- function(x){
  Mean = round(mean(x, na.rm=T),2)
  Sd = round(sd(x, na.rm=T),2)
  return(paste0(Mean, '+-',Sd))
}
pure_df <- function(x){
  out = x[!is.na(x)]
  return(out)
}

fiber_vertex_file = 'BrainGraphRDS/CFS_fiber_weighted_subject_vertex.csv'
fiber_vertex = fread(fiber_vertex_file)
fiber_vertex = fiber_vertex[threshold==0.25, ]
fiber_vertex_146 = fiber_vertex[participant_id %in% clinic_info_73$participant_id, ]

roi_index = 198
aa_p = fiber_vertex_146[region == brainnetome_surface$name[roi_index]&all_group=='patient_before', .(degree)]
Mean_Sd(aa_p[[1]])
aa_h = fiber_vertex_146[region == brainnetome_surface$name[roi_index]&all_group=='health_before', .(degree)]
Mean_Sd(aa_h[[1]])
t.test(aa_p[[1]], aa_h[[1]], var.equal = T)

roi_index = 115
aa_p = fiber_vertex_146[region == brainnetome_surface$name[roi_index]&all_group=='patient_before', .(dist)]
Mean_Sd(aa_p[[1]])
aa_h = fiber_vertex_146[region == brainnetome_surface$name[roi_index]&all_group=='health_before', .(dist)]
Mean_Sd(aa_h[[1]])
t.test(pure_df(aa_p[[1]]), pure_df(aa_h[[1]]), var.equal = T)

roi_index = 29
aa_p = fiber_vertex_146[region == brainnetome_surface$name[roi_index]&all_group=='patient_before', .(E.nodal.wt)]
Mean_Sd(aa_p[[1]])
aa_h = fiber_vertex_146[region == brainnetome_surface$name[roi_index]&all_group=='health_before', .(E.nodal.wt)]
Mean_Sd(aa_h[[1]])
t.test(aa_p[[1]], aa_h[[1]], var.equal = T)

##############################################
library(ggplot2)
indicator = 'degree'
roi_index = 198
roi_name = brainnetome_surface$name[roi_index]
group_name = c('patient_before', 'health_before')
select_ind = c(indicator, 'participant_id', 'all_group')
select_fiber = fiber_vertex_146[all_group%in%group_name&region==roi_name, ..select_ind]
# refer to DBS group
# two groups: blue #8cb5c1, red f1afa0
select_fiber$all_group <- factor(select_fiber$all_group, levels = c('patient_before', 'health_before'))
ggplot(select_fiber, aes(x=all_group, y=degree, fill=all_group)) +
  geom_violin(trim=F, color='#406AF2', size=1.5, alpha=.5) +
  geom_boxplot(width=0.05, fill='white') + 
  theme_classic() +
  scale_y_continuous(expand = c(0,0), limits = c(20, 100)) +
  scale_x_discrete(labels = c("health_before" = "HC", "patient_before" = "CFS")) +
  labs(x='', y='Nodal Degree') +
  theme(
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    text = element_text(family = "Times New Roman"),
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    axis.text = element_text(size=24),
    legend.position = 'none')

##############################################
# distance
indicator = 'dist'
roi_index = 115
roi_name = brainnetome_surface$name[roi_index]
group_name = c('patient_before', 'health_before')
select_ind = c(indicator, 'participant_id', 'all_group')
select_fiber = fiber_vertex_146[all_group%in%group_name&region==roi_name, ..select_ind]
# refer to DBS group
# two groups: blue #8cb5c1, red f1afa0
select_fiber$all_group <- factor(select_fiber$all_group, levels = c('patient_before', 'health_before'))
ggplot(select_fiber, aes(x=all_group, y=dist, fill=all_group)) +
  geom_violin(trim=F, color='#406AF2', size=1.5, alpha=.5) +
  geom_boxplot(width=0.05, fill='white') + 
  theme_classic() +
  scale_y_continuous(expand = c(0,0), limits = c(10, 80)) +
  scale_x_discrete(labels = c("health_before" = "HC", "patient_before" = "CFS")) +
  labs(x='', y='Nodal Distance') +
  theme(
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    text = element_text(family = "Times New Roman"),
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    axis.text = element_text(size=24),
    legend.position = 'none')


##############################################
# nodal efficiency
indicator = 'E.nodal.wt'
roi_index = 29
roi_name = brainnetome_surface$name[roi_index]
group_name = c('patient_before', 'health_before')
select_ind = c(indicator, 'participant_id', 'all_group')
select_fiber = fiber_vertex_146[all_group%in%group_name&region==roi_name, ..select_ind]
# refer to DBS group
# two groups: blue #8cb5c1, red f1afa0
select_fiber$all_group <- factor(select_fiber$all_group, levels = c('patient_before', 'health_before'))
select_fiber = select_fiber[E.nodal.wt>1, ]
ggplot(select_fiber, aes(x=all_group, y=E.nodal.wt, fill=all_group)) +
  geom_violin(trim=F, color='#406AF2', size=1.5, alpha=.5) +
  geom_boxplot(width=0.05, fill='white') + 
  theme_classic() +
  scale_y_continuous(expand = c(0,0), limits = c(1, 3)) +
  scale_x_discrete(labels = c("health_before" = "HC", "patient_before" = "CFS")) +
  labs(x='', y='Weighted Nodal Efficiency') +
  theme(
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    text = element_text(family = "Times New Roman"),
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    axis.text = element_text(size=24),
    legend.position = 'none')

##############################################
#### pre- & post-
# nodal degree
clinic_info_73 = fread('subs_146_info.csv')
load('brainnetome_surface.rda')
Mean_Sd <- function(x){
  Mean = round(mean(x, na.rm=T),2)
  Sd = round(sd(x, na.rm=T),2)
  return(paste0(Mean, '+-',Sd))
}

fiber_vertex_file = 'BrainGraphRDS/CFS_fiber_weighted_subject_vertex.csv'
fiber_vertex = fread(fiber_vertex_file)
fiber_vertex = fiber_vertex[threshold==0.25, ]
fiber_vertex_146 = fiber_vertex[participant_id %in% clinic_info_73$participant_id, ]

compute_diff <- function(df, group_type, indicator){ # df: datatable, group:patient, health
  group_df = df[group==group_type, ]
  setkey(group_df, all_group, subject)
  sub_len = length(group_df$subject)
  sub_data = group_df[[indicator]]
  return(sub_data[1:(sub_len/2)] - sub_data[(sub_len/2+1):sub_len])
}

roi_index = 198
roi_name = brainnetome_surface$name[roi_index]
aa_p = fiber_vertex_146[region == roi_name&all_group=='patient_after', .(degree)]
Mean_Sd(aa_p[[1]])
aa_h = fiber_vertex_146[region == roi_name&all_group=='health_after', .(degree)]
Mean_Sd(aa_h[[1]])
roi_df <- merge( clinic_info_73, 
              fiber_vertex_146[region == roi_name, .(degree, participant_id)], 
              by='participant_id' )
roi_diff <- compute_diff(roi_df, 'health', 'degree')
Mean_Sd(roi_diff)
wilcox.test(roi_diff)


roi_index = 115
roi_name = brainnetome_surface$name[roi_index]
aa_p = fiber_vertex_146[region == roi_name&all_group=='patient_after', .(dist)]
Mean_Sd(aa_p[[1]])
aa_h = fiber_vertex_146[region == roi_name&all_group=='health_after', .(dist)]
Mean_Sd(aa_h[[1]])
roi_df <- merge( clinic_info_73, 
                 fiber_vertex_146[region == roi_name, .(dist, participant_id)], 
                 by='participant_id' )
roi_diff <- compute_diff(roi_df, 'patient', 'dist')
Mean_Sd(roi_diff)
wilcox.test(roi_diff)
roi_df <- merge( clinic_info_73, 
                 fiber_vertex_146[region == roi_name, .(dist, participant_id)], 
                 by='participant_id' )
roi_diff <- compute_diff(roi_df, 'health', 'dist')
Mean_Sd(roi_diff)
wilcox.test(roi_diff)


roi_index = 29
roi_name = brainnetome_surface$name[roi_index]
aa_p = fiber_vertex_146[region == roi_name&all_group=='patient_after', .(E.nodal.wt)]
Mean_Sd(aa_p[[1]])
aa_h = fiber_vertex_146[region == roi_name&all_group=='health_after', .(E.nodal.wt)]
Mean_Sd(aa_h[[1]])
roi_df <- merge( clinic_info_73, 
                 fiber_vertex_146[region == roi_name, .(E.nodal.wt, participant_id)], 
                 by='participant_id' )
roi_diff <- compute_diff(roi_df, 'patient', 'E.nodal.wt')
Mean_Sd(roi_diff)
t.test(roi_diff)
roi_df <- merge( clinic_info_73, 
                 fiber_vertex_146[region == roi_name, .(E.nodal.wt, participant_id)], 
                 by='participant_id' )
roi_diff <- compute_diff(roi_df, 'health', 'E.nodal.wt')
Mean_Sd(roi_diff)
t.test(roi_diff)

###### correlation
diff_correlation <- function(df, group_type, indicator, method){ # df: datatable, group:patient, health
  group_df = df[group==group_type, ]
  setkey(group_df, all_group, subject)
  sub_len = length(group_df$subject)
  sub_net = group_df[[indicator]]
  sub_net_diff = sub_net[1:(sub_len/2)] - sub_net[(sub_len/2+1):sub_len]
  sub_fs14 = group_df[['FS14']]
  sub_fs14_diff = sub_fs14[1:(sub_len/2)] - sub_fs14[(sub_len/2+1):sub_len]
  sub_fs14_cor = cor.test(sub_net_diff, sub_fs14_diff, method = method)
  return(cat(paste('\n',method,
                   '    \nStatistic value',round(sub_fs14_cor$statistic,2),
                   '    \nR value',round(sub_fs14_cor$estimate,2),
                   '    \nP value',round(sub_fs14_cor$p.value,3))))
}

roi_index = 198
roi_name = brainnetome_surface$name[roi_index]
roi_df <- merge( clinic_info_73, 
                 fiber_vertex_146[region == roi_name, .(degree, participant_id)], 
                 by='participant_id' )
diff_correlation(roi_df, 'patient', 'degree', 'spearman')
diff_correlation(roi_df, 'health', 'degree', 'spearman')

roi_index = 115
roi_name = brainnetome_surface$name[roi_index]
roi_df <- merge( clinic_info_73, 
                 fiber_vertex_146[region == roi_name, .(dist, participant_id)], 
                 by='participant_id' )
diff_correlation(roi_df, 'patient', 'dist', 'spearman')
diff_correlation(roi_df, 'health', 'dist', 'spearman')

roi_index = 29
roi_name = brainnetome_surface$name[roi_index]
roi_df <- merge( clinic_info_73, 
                 fiber_vertex_146[region == roi_name, .(E.nodal.wt, participant_id)], 
                 by='participant_id' )
diff_correlation(roi_df, 'patient', 'E.nodal.wt', 'pearson')
diff_correlation(roi_df, 'health', 'E.nodal.wt', 'pearson')

#############################################################
## all brain

#### each region
## degree
diff_correlation_df <- function(df, group_type, indicator){ 
  group_df = df[group==group_type, ]
  setkey(group_df, all_group, subject)
  sub_len = length(group_df$subject)
  sub_net = group_df[[indicator]]
  sub_net_diff = sub_net[1:(sub_len/2)] - sub_net[(sub_len/2+1):sub_len]
  sub_fs14 = group_df[['FS14']]
  sub_fs14_diff = sub_fs14[1:(sub_len/2)] - sub_fs14[(sub_len/2+1):sub_len]
  sub_out = list()
  sub_out[[1]] = sub_net_diff
  sub_out[[2]] = sub_fs14_diff
  return(sub_out)
}

indicator = 'E.nodal.wt'
roi_index = 29
roi_name = brainnetome_surface$name[roi_index]
select_col = c(indicator, 'participant_id')
roi_df <- merge( clinic_info_73, fiber_vertex_146[region == roi_name, ..select_col], 
                 by='participant_id' )
diff_list <- diff_correlation_df(roi_df, 'health', indicator)
diff_df = data.table(net=diff_list[[1]], fs14=diff_list[[2]])

ggplot(diff_df, aes(x=fs14, y=net)) +
  geom_point(size=4) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth=4) +
  theme_test() + 
  labs(x='FS14 Difference', y='Weighted Nodal Efficiency', 
       title='HC, r = 0.41, P = 0.011') +
  theme(
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    text = element_text(family = "Times New Roman"),
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    axis.text = element_text(size = 24),
    plot.title = element_text(size = 32, hjust = 0.1))  




