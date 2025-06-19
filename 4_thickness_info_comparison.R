library(data.table)
options(warn = -1) 
library(coin)
############################## vertex
clinic_info_73 = fread('subs_146_info.csv')
load('brainnetome_surface.rda')

func_vertex_file = 'BrainGraphRDS/CFS_coupling_residual_scn_vertex.csv'
func_vertex = fread(func_vertex_file)

#[1] "density"       "region"        "lobe"          "hemi"          "gyrus"        
#[6] "Yeo_7network"  "Yeo_17network" "degree"        "comp"          "asymm"        
#[11] "dist"          "dist.strength" "knn"           "Lp"            "btwn.cent"    
#[16] "hubs"          "ev.cent"       "lev.cent"      "k.core"        "transitivity" 
#[21] "E.local"       "E.nodal"       "vulnerability" "eccentricity"  "comm"         
#[26] "GC"            "PC"            "z.score"       "atlas"         "modality"     
#[31] "all_group"     "threshold"  

select_col = 'dist'
net_df = func_vertex
corr_between_net_behavior_vertex <- function(net_df, select_col){
  net_select_col = c(select_col, 'region', 'all_group', 'threshold' )
  net_select_df = net_df[, ..net_select_col]
  colnames(net_select_df) <- c('net_value', 'region', 'all_group', 'threshold')
  merge_df = merge(brainnetome_surface[, .(name, index)], net_select_df, by.x = 'name', by.y = 'region')    
  setkey(merge_df, 'index')  
  patient_before = merge_df[all_group=='patient_before', ]  
  patient_after = merge_df[all_group=='patient_after', ]  
  health_before = merge_df[all_group=='health_before', ]  
  health_after = merge_df[all_group=='health_after', ]  
  
  # all regions
  all_regions = brainnetome_surface$name
  p_before_h_before_z = vector()
  p_before_h_before_p = vector()
  p_before_p_after_z = vector()
  p_before_p_after_p = vector()
  h_before_h_after_z = vector()
  h_before_h_after_p = vector()

  for (i in seq_along(all_regions)){
    i_name = all_regions[i]
    ##### p_before, h_before
    df = rbind(patient_before[name==i_name, ], health_before[name==i_name, ])
    df$all_group = as.factor(df$all_group)
    region_perm <- independence_test(net_value ~ all_group, data = df, distribution = approximate(B = 1000))
    p_before_h_before_p[i] = pvalue(region_perm)
    p_before_h_before_z[i] = statistic(region_perm)
    ##### p_before, p_after
    df = rbind(patient_before[name==i_name, ], patient_after[name==i_name, ])
    df$all_group = as.factor(df$all_group)
    region_perm <- independence_test(net_value ~ all_group, data = df, distribution = approximate(B = 1000))
    p_before_p_after_p[i] = pvalue(region_perm)
    p_before_p_after_z[i] = statistic(region_perm) 
    ##### h_before, h_after
    df = rbind(health_before[name==i_name, ], health_after[name==i_name, ])
    df$all_group = as.factor(df$all_group)
    region_perm <- independence_test(net_value ~ all_group, data = df, distribution = approximate(B = 1000))
    h_before_h_after_p[i] = pvalue(region_perm)
    h_before_h_after_z[i] = statistic(region_perm)     
  }

  out_df = data.table(
    p_before_h_before_p = p_before_h_before_p,
    p_before_p_after_p = p_before_p_after_p,
    h_before_h_after_p = h_before_h_after_p
  )
  return(out_df)
}
my_df = corr_between_net_behavior_vertex(net_df, select_col)
fwrite(my_df, file='dist.csv')

##############################################
library(ciftiTools)
ciftiTools.setOption('wb_path', '/usr/local/workbench')
file = '/home/clancy/TemplateFlow/BN_Atlas_freesurfer/100307/fsaverage_LR32k/100307.BN_Atlas.32k_fs_LR.dlabel.nii'
bn = read_cifti(file)
bn_lh = bn$data$cortex_left[,1]
bn_rh = bn$data$cortex_right[,1]
bn_rh_new = bn_rh - 210
bn_rh_new[bn_rh_new<0] = 0
bn_both = c(bn_lh, bn_rh_new)
generate_xifti <- function(empty_num, input_data){
  for (i in seq(210)){
    empty_num[bn_both==i] <- input_data[i]
  }
  label_xifti <- newdata_xifti(bn, empty_num)
  return(label_xifti)
}
##############################################
##########################
load('brainnetome_surface.rda')
## thickness
func_vertex_file = 'BrainGraphRDS/CFS_thickness_scn_vertex.csv'
func_vertex = fread(func_vertex_file)
func_vertex = func_vertex[threshold==0.25&all_group=='patient_after', ]
net_indicator = 'degree'
df_average = merge(func_vertex, brainnetome_surface[, .(name, index)], 
                   by.x = 'region', by.y = 'name')
setkey(df_average, index)
empty_num = rep(NA, length(bn_both))
input_data = df_average[[net_indicator]] # 210 values
label_xifti = generate_xifti(empty_num, input_data)

library(RColorBrewer)
view_xifti_surface(label_xifti,  zlim=c(30,80), color_mode='diverging',colors=(brewer.pal(10, 'RdPu')) )


view_xifti_surface(label_xifti, zlim=c(-7, 7), color_mode='diverging', 
                   colors=(brewer.pal(10, 'RdPu')),
                   fname='Degree_health_diff.png')

Mean_Sd <- function(x){
  Mean = round(mean(x),2)
  Sd = round(sd(x),2)
  return(paste0(Mean, '+-',Sd))
}

#############################################################
net_indicator = 'degree'
roi_index = 198
group_name = 'health_after'
func_vertex_file = 'BrainGraphRDS/CFS_coupling_residual_scn_vertex.csv'
func_vertex = fread(func_vertex_file)
func_vertex = func_vertex[threshold==0.25, ]
roi_name = brainnetome_surface[index==roi_index, .(name)][[1]]
select_col = c(net_indicator, 'all_group', 'threshold')
func_vertex[region==roi_name&all_group==group_name, ..select_col ]


df = rbind(patient_before[name==i_name, ], health_before[name==i_name, ])
df$all_group = as.factor(df$all_group)
region_perm <- independence_test(net_value ~ all_group, data = df, distribution = approximate(B = 1000))



indicator = 'degree'
select_ind = c(indicator, 'participant_id', 'region')
select_func = func_vertex_146[all_group=='health_before', ..select_ind]
df_average = merge(brainnetome_surface, select_func[, .(Value=mean(degree)), by=.(region)], 
                   by.x = 'name', by.y = 'region')
# patient_before: max 74.27778, health_before: max 71.10811
setkey(df_average, index)
empty_num = rep(NA, length(bn_both))
input_data = df_average$Value # 210 values
label_xifti = generate_xifti(empty_num, input_data)
view_xifti_surface(label_xifti, zlim=c(30,70), color_mode='diverging',
                   fname = 'Degree_health_before.png')
#bn_193 <- transform_xifti(bn, function(x){ifelse(x==193, 193, NA)})
#plot(bn_193, borders='black', title='BNA 193', fname='bn_193.png')
library(ggplot2)
indicator = 'degree'
roi_index = 193
roi_name = brainnetome_surface$name[roi_index]
group_name = c('patient_before', 'health_before')
select_ind = c(indicator, 'participant_id', 'all_group')
select_func = func_vertex_146[all_group%in%group_name&region==roi_name, ..select_ind]
# refer to DBS group
# two groups: blue #8cb5c1, red f1afa0
select_func$all_group <- factor(select_func$all_group, levels = c('patient_before', 'health_before'))
ggplot(select_func, aes(x=all_group, y=degree, fill=all_group)) +
  geom_violin(trim=F, color='#406AF2', size=1.5, alpha=.5) +
  geom_boxplot(width=0.05, fill='white') + 
  theme_classic() +
  scale_y_continuous(expand = c(0,0), limits = c(0, 80)) +
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
select_ind = c(indicator, 'participant_id', 'region')
select_func = func_vertex_146[all_group=='patient_before', ..select_ind]
df_average = merge(brainnetome_surface, select_func[, .(Value=mean(dist)), by=.(region)], 
                   by.x = 'name', by.y = 'region')
# patient_before: max 91.60951, health_before: max 91.60968
setkey(df_average, index)
empty_num = rep(NA, length(bn_both))
input_data = df_average$Value # 210 values
label_xifti = generate_xifti(empty_num, input_data)
#view_xifti_surface(label_xifti, zlim=c(50,80), color_mode='diverging')
view_xifti_surface(label_xifti, zlim=c(50,80), color_mode='diverging',
                   fname = 'Dist_patient_before.png')
#bn_192 <- transform_xifti(bn, function(x){ifelse(x==193 | x==402, 193, NA)}) # right add 210
#plot(bn_192, borders='black', title='BNA 192 193', fname='bn_192_193.png')

indicator = 'dist'
roi_index = 193
roi_name = brainnetome_surface$name[roi_index]
group_name = c('patient_before', 'health_before')
select_ind = c(indicator, 'participant_id', 'all_group')
select_func = func_vertex_146[all_group%in%group_name&region==roi_name, ..select_ind]
select_func$all_group <- factor(select_func$all_group, levels = c('patient_before', 'health_before'))
ggplot(select_func, aes(x=all_group, y=dist, fill=all_group)) +
  geom_violin(trim=F, color='#406AF2', size=1.5, alpha=.5) +
  geom_boxplot(width=0.05, fill='white') + 
  theme_classic() +
  scale_y_continuous(expand = c(0,0), limits = c(10, 100)) +
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
indicator = 'E.nodal'
select_ind = c(indicator, 'participant_id', 'region')
select_func = func_vertex_146[all_group=='health_before', ..select_ind]
df_average = merge(brainnetome_surface, select_func[, .(Value=mean(E.nodal)), by=.(region)], 
                   by.x = 'name', by.y = 'region')
# patient_before: max 0.6670875, health_before: max 0.6602116
setkey(df_average, index)
empty_num = rep(NA, length(bn_both))
input_data = df_average$Value # 210 values
label_xifti = generate_xifti(empty_num, input_data)
view_xifti_surface(label_xifti, zlim=c(0.55,0.65), color_mode='diverging',
                   fname = 'Enodal_health_before.png')

view_xifti_surface(label_xifti, zlim=c(0,80), color_mode='diverging',
                   fname = 'Dist_patient_before.png')
bn_192 <- transform_xifti(bn, function(x){ifelse(x==189 | x==400 | x==193, x, NA)})
plot(bn_192, borders=FALSE, title='BNA 189 190 193', 
     color_mode='qualitative', fname='bn_189_190_193.png')

indicator = 'E.nodal'
roi_index = 193
roi_name = brainnetome_surface$name[roi_index]
group_name = c('patient_before', 'health_before')
select_ind = c(indicator, 'participant_id', 'all_group')
select_func = func_vertex_146[all_group%in%group_name&region==roi_name, ..select_ind]
select_func$all_group <- factor(select_func$all_group, levels = c('patient_before', 'health_before'))
ggplot(select_func, aes(x=all_group, y=E.nodal, fill=all_group)) +
  geom_violin(trim=F, color='#406AF2', size=1.5, alpha=.5) +
  geom_boxplot(width=0.05, fill='white') + 
  theme_classic() +
  scale_y_continuous(expand = c(0,0), limits = c(0.35, 0.75)) +
  scale_x_discrete(labels = c("health_before" = "HC", "patient_before" = "CFS")) +
  labs(x='', y='Nodal Efficiency') +
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
  Mean = round(mean(x),2)
  Sd = round(sd(x),2)
  return(paste0(Mean, '+-',Sd))
}

func_vertex_file = 'BrainGraphRDS/CFS_func_weighted_subject_vertex.csv'
func_vertex = fread(func_vertex_file)
func_vertex = func_vertex[threshold==0.25, ]
func_vertex_146 = func_vertex[participant_id %in% clinic_info_73$participant_id, ]

compute_diff <- function(df, group_type, indicator){ # df: datatable, group:patient, health
  group_df = df[group==group_type, ]
  setkey(group_df, all_group, subject)
  sub_len = length(group_df$subject)
  sub_data = group_df[[indicator]]
  return(sub_data[1:(sub_len/2)] - sub_data[(sub_len/2+1):sub_len])
}

roi_index = 193
roi_name = brainnetome_surface$name[roi_index]
aa_p = func_vertex_146[region == roi_name&all_group=='patient_after', .(degree)]
Mean_Sd(aa_p[[1]])
aa_h = func_vertex_146[region == roi_name&all_group=='health_after', .(degree)]
Mean_Sd(aa_h[[1]])
roi_df <- merge( clinic_info_73, 
                 func_vertex_146[region == roi_name, .(degree, participant_id)], 
                 by='participant_id' )
roi_diff <- compute_diff(roi_df, 'health', 'degree')
Mean_Sd(roi_diff)
wilcox.test(roi_diff)


roi_index = 193
roi_name = brainnetome_surface$name[roi_index]
aa_p = func_vertex_146[region == roi_name&all_group=='patient_after', .(dist)]
Mean_Sd(aa_p[[1]])
aa_h = func_vertex_146[region == roi_name&all_group=='health_after', .(dist)]
Mean_Sd(aa_h[[1]])
roi_df <- merge( clinic_info_73, 
                 func_vertex_146[region == roi_name, .(dist, participant_id)], 
                 by='participant_id' )
roi_diff <- compute_diff(roi_df, 'patient', 'dist')
Mean_Sd(roi_diff)
wilcox.test(roi_diff)
roi_df <- merge( clinic_info_73, 
                 func_vertex_146[region == roi_name, .(dist, participant_id)], 
                 by='participant_id' )
roi_diff <- compute_diff(roi_df, 'health', 'dist')
Mean_Sd(roi_diff)
wilcox.test(roi_diff)


roi_index = 189
roi_name = brainnetome_surface$name[roi_index]
aa_p = func_vertex_146[region == roi_name&all_group=='patient_after', .(E.nodal)]
Mean_Sd(aa_p[[1]])
aa_h = func_vertex_146[region == roi_name&all_group=='health_after', .(E.nodal)]
Mean_Sd(aa_h[[1]])
roi_df <- merge( clinic_info_73, 
                 func_vertex_146[region == roi_name, .(E.nodal, participant_id)], 
                 by='participant_id' )
roi_diff <- compute_diff(roi_df, 'patient', 'E.nodal')
Mean_Sd(roi_diff)
t.test(roi_diff)
roi_df <- merge( clinic_info_73, 
                 func_vertex_146[region == roi_name, .(E.nodal, participant_id)], 
                 by='participant_id' )
roi_diff <- compute_diff(roi_df, 'health', 'E.nodal')
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

roi_index = 193
roi_name = brainnetome_surface$name[roi_index]
roi_df <- merge( clinic_info_73, 
                 func_vertex_146[region == roi_name, .(degree, participant_id)], 
                 by='participant_id' )
diff_correlation(roi_df, 'patient', 'degree', 'spearman')
diff_correlation(roi_df, 'health', 'degree', 'spearman')

roi_index = 193
roi_name = brainnetome_surface$name[roi_index]
roi_df <- merge( clinic_info_73, 
                 func_vertex_146[region == roi_name, .(E.nodal, participant_id)], 
                 by='participant_id' )
diff_correlation(roi_df, 'patient', 'E.nodal', 'pearson')
diff_correlation(roi_df, 'health', 'E.nodal', 'pearson')

#############################################################
## all brain
indicator = 'degree'
p_before_id = clinic_info_73[all_group=='patient_before', .(participant_id, subject)]
p_after_id = clinic_info_73[all_group=='patient_after', .(participant_id, subject)]
setkey(p_before_id, subject, participant_id)
setkey(p_after_id, subject, participant_id)
subj_list = matrix(NA, nrow=length(p_before_id$subject), ncol=210)
for (i in seq_along(p_before_id$subject)){
  i_subject = p_before_id$subject[i]
  i_id_before = p_before_id[subject==i_subject,.(participant_id)][[1]]
  i_id_after = p_after_id[subject==i_subject,.(participant_id)][[1]]
  selec_col = c('participant_id', 'region', indicator)
  i_after_df = merge(brainnetome_surface, func_vertex_146[participant_id==i_id_after, ..selec_col],
                     by.x = 'name', by.y = 'region')
  setkey(i_after_df, index)
  i_before_df = merge(brainnetome_surface, func_vertex_146[participant_id==i_id_before, ..selec_col],
                      by.x = 'name', by.y = 'region')
  setkey(i_before_df, index)
  i_diff = i_after_df[[indicator]] - i_before_df[[indicator]]
  subj_list[i,] = i_diff
}

df_average <- colMeans(subj_list)
# patient: -10.94444 ~ 9.416667, health_before: -7.027027 ~ 6.972973
empty_num = rep(NA, length(bn_both))
input_data = df_average # 210 values
label_xifti = generate_xifti(empty_num, input_data)
library(RColorBrewer)
view_xifti_surface(label_xifti, zlim=c(-7, 7), color_mode='diverging', 
                   colors=rev(brewer.pal(10, 'RdBu')),
                   fname='Degree_health_diff.png')

## dist
indicator = 'dist'
p_before_id = clinic_info_73[all_group=='health_before', .(participant_id, subject)]
p_after_id = clinic_info_73[all_group=='health_after', .(participant_id, subject)]
setkey(p_before_id, subject, participant_id)
setkey(p_after_id, subject, participant_id)
subj_list = matrix(NA, nrow=length(p_before_id$subject), ncol=210)
for (i in seq_along(p_before_id$subject)){
  i_subject = p_before_id$subject[i]
  i_id_before = p_before_id[subject==i_subject,.(participant_id)][[1]]
  i_id_after = p_after_id[subject==i_subject,.(participant_id)][[1]]
  selec_col = c('participant_id', 'region', indicator)
  i_after_df = merge(brainnetome_surface, func_vertex_146[participant_id==i_id_after, ..selec_col],
                     by.x = 'name', by.y = 'region')
  setkey(i_after_df, index)
  i_before_df = merge(brainnetome_surface, func_vertex_146[participant_id==i_id_before, ..selec_col],
                      by.x = 'name', by.y = 'region')
  setkey(i_before_df, index)
  i_diff = i_after_df[[indicator]] - i_before_df[[indicator]]
  subj_list[i,] = i_diff
}

df_average <- colMeans(subj_list)
# patient: -2.395498 ~ 10.00974, health_before: -3.928684 ~ 5.388388
empty_num = rep(NA, length(bn_both))
input_data = df_average # 210 values
label_xifti = generate_xifti(empty_num, input_data)
library(RColorBrewer)
view_xifti_surface(label_xifti, zlim=c(-7, 7), color_mode='diverging', 
                   colors=rev(brewer.pal(10, 'RdBu')),
                   fname='Dist_health_diff.png')

## E.nodal
indicator = 'E.nodal'
p_before_id = clinic_info_73[all_group=='health_before', .(participant_id, subject)]
p_after_id = clinic_info_73[all_group=='health_after', .(participant_id, subject)]
setkey(p_before_id, subject, participant_id)
setkey(p_after_id, subject, participant_id)
subj_list = matrix(NA, nrow=length(p_before_id$subject), ncol=210)
for (i in seq_along(p_before_id$subject)){
  i_subject = p_before_id$subject[i]
  i_id_before = p_before_id[subject==i_subject,.(participant_id)][[1]]
  i_id_after = p_after_id[subject==i_subject,.(participant_id)][[1]]
  selec_col = c('participant_id', 'region', indicator)
  i_after_df = merge(brainnetome_surface, func_vertex_146[participant_id==i_id_after, ..selec_col],
                     by.x = 'name', by.y = 'region')
  setkey(i_after_df, index)
  i_before_df = merge(brainnetome_surface, func_vertex_146[participant_id==i_id_before, ..selec_col],
                      by.x = 'name', by.y = 'region')
  setkey(i_before_df, index)
  i_diff = i_after_df[[indicator]] - i_before_df[[indicator]]
  subj_list[i,] = i_diff
}

df_average <- colMeans(subj_list)
# patient: -0.02819865 ~ 0.02177477, health_before: -0.01163843 ~ 0.01976378
empty_num = rep(NA, length(bn_both))
input_data = df_average # 210 values
label_xifti = generate_xifti(empty_num, input_data)
library(RColorBrewer)
view_xifti_surface(label_xifti, zlim=c(-0.02, 0.020), color_mode='diverging', 
                   colors=rev(brewer.pal(10, 'RdBu')),
                   fname='Enodal_health_diff.png')

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

indicator = 'E.nodal'
roi_index = 193
roi_name = brainnetome_surface$name[roi_index]
select_col = c(indicator, 'participant_id')
roi_df <- merge( clinic_info_73, func_vertex_146[region == roi_name, ..select_col], 
                 by='participant_id' )
diff_list <- diff_correlation_df(roi_df, 'patient', indicator)
diff_df = data.table(net=diff_list[[1]], fs14=diff_list[[2]])

ggplot(diff_df, aes(x=fs14, y=net)) +
  geom_point(size=4) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth=4) +
  theme_test() + 
  labs(x='FS14 Difference', y='Nodal Efficiency Difference', 
       title='CFS, r = 0.38, P = 0.021') +
  theme(
    plot.margin = unit(c(1, 1, 1, 1), "cm"),
    text = element_text(family = "Times New Roman"),
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    axis.text = element_text(size = 24),
    plot.title = element_text(size = 32, hjust = 0.1))  



