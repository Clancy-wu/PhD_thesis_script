library(data.table)
clinic_info_file = 'subs_146_info.csv'
whole_brain_file = 'Measurements/whole_brain.csv'
alff_file = 'Measurements/alff_results.csv'
thickness_file = 'Measurements/CFS_thickness_mat_exclude_Subject_019.csv'
coupling_file = 'Measurements/CFS_coupling_mat_exclude_Subject_019.csv'
# "ROI-29" "ROI-115" "ROI-189" "ROI-190" "ROI-192" "ROI-193" "ROI-198"  
Mean_Sd <- function(x){
  Mean = round(mean(x),2)
  Sd = round(sd(x),2)
  return(paste0(Mean, '+-',Sd))
}
####################################

# 1. alff
thickness = fread(thickness_file)
thickness_sum = rowSums(thickness[,2:211])
thickness$Sum = thickness_sum
thickness_col = c('subject', "ROI-29", "ROI-115", "ROI-189", "ROI-190", "ROI-192", "ROI-193", "ROI-198", 'Sum')
thickness_df = thickness[, ..thickness_col]
clinic_info = fread(clinic_info_file)
alff_results = merge(clinic_info, thickness_df, by.x = 'participant_id', by.y = 'subject')
#alff_results = fread(alff_file)
roi_indexs = c("ROI-29", "ROI-115", "ROI-189", "ROI-190", "ROI-192", "ROI-193", "ROI-198")
df = alff_results
roi_index = roi_indexs[4]
# df: participant_id, subject, all_group, roi_index
#### whole brain
#whole_brain = fread(whole_brain_file)
#df = merge(df, whole_brain, by='participant_id')
df = df[subject != 'Subject_019', ]
#select_col = c('participant_id', 'subject', 'group', 'treatment', 'all_group', 'brain_size', roi_index, 'FS14' )
select_col = c('participant_id', 'subject', 'group', 'treatment', 'all_group', 'Sum', roi_index, 'FS14' )
df_prepare = df[, ..select_col]
colnames(df_prepare) <- c('participant_id', 'subject', 'group', 'treatment', 'all_group', 'brain_size', 'Value', 'FS14' )
## 1. baseline diff
df_baseline = df_prepare[treatment=='before', ]
tapply(df_baseline$Value, df_baseline$group, Mean_Sd)
model = aov(Value ~ group + brain_size, data=df_baseline)
summary(model)
#####################################################
df_patient = df_prepare[group=='patient', ]
library(afex)
library(effectsize)
mod <- afex::lmer_alt(Value ~ treatment + 1 + (1|subject), data=df_patient, method='PB')
summary(mod)


#####################################################
# 2. patient pre- post-
p_before = df_prepare[all_group=='health_before', ]
p_after = df_prepare[all_group=='health_after', ]
setkey(p_before, subject); setkey(p_after, subject)
Mean_Sd(p_after$Value)
p_diff = p_after$Value - p_before$Value
Mean_Sd(p_diff)
t.test(p_diff)
cor.test(p_after$Value - p_before$Value, p_after$FS14 - p_before$FS14, method = 'pearson')








## 3. health pre- post-
h_before = df_prepare[all_group=='health_before', ]
h_after = df_prepare[all_group=='health_after', ]
setkey(h_before, subject); setkey(h_after, subject)
wilcox.test(h_after$Value - h_before$Value)
cor.test(h_after$Value - h_before$Value, h_after$FS14 - h_before$FS14, method = 'pearson') 





