library(data.table)
#### clinical info
clinic_info = fread('all_subs_152_info.csv')
clinic_before = clinic_info[treatment=='before', ]
Mean_Sd <- function(x){
  Mean = round(mean(x),2)
  Sd = round(sd(x),2)
  return(paste0(Mean, '+-',Sd))
}
Mean_Sd(clinic_before[group=='patient', (age)])
Mean_Sd(clinic_before[group=='health', (age)])
tapply(clinic_before$age, clinic_before$group, shapiro.test)
wilcox.test(clinic_before$age ~ clinic_before$group)

table(clinic_before[group=='patient', (sex)])
table(clinic_before[group=='health', (sex)])
prop.test(c(11,12), c(39,38))

Mean_Sd(clinic_before[group=='patient', (BMI)])
Mean_Sd(clinic_before[group=='health', (BMI)])
tapply(clinic_before$BMI, clinic_before$group, shapiro.test)
t.test(clinic_before$BMI ~ clinic_before$group, var.equal=T)

Mean_Sd(clinic_before[group=='patient', (disease_duration)])

Mean_Sd(clinic_before[group=='patient', (FS14)])
Mean_Sd(clinic_before[group=='health', (FS14)])
tapply(clinic_before$FS14, clinic_before$group, shapiro.test)
wilcox.test(clinic_before$FS14 ~ clinic_before$group)

Mean_Sd(clinic_before[group=='patient', (SF36)])
Mean_Sd(clinic_before[group=='health', (SF36)])
tapply(clinic_before$SF36, clinic_before$group, shapiro.test)
wilcox.test(clinic_before$SF36 ~ clinic_before$group)

Mean_Sd(clinic_before[group=='patient', (PSQI)])
Mean_Sd(clinic_before[group=='health', (PSQI)])
tapply(clinic_before$PSQI, clinic_before$group, shapiro.test)
t.test(clinic_before$PSQI ~ clinic_before$group, var.equal=T)

Mean_Sd(clinic_before[group=='patient', (HAMD)])
Mean_Sd(clinic_before[group=='health', (HAMD)])
tapply(clinic_before$HAMD, clinic_before$group, shapiro.test)
wilcox.test(clinic_before$HAMD ~ clinic_before$group)

Mean_Sd(clinic_before[group=='patient', (PRI)])
Mean_Sd(clinic_before[group=='health', (PRI)])
tapply(clinic_before$PRI, clinic_before$group, shapiro.test)
wilcox.test(clinic_before$PRI ~ clinic_before$group)

#### 
Mean_Sd(clinic_info[treatment=='after' & group=='patient', (PRI)])
Mean_Sd(clinic_info[treatment=='after' & group=='health', (PRI)])

## LME
library(afex)
library(effectsize)
mod <- afex::lmer_alt(PRI ~ group + treatment + group*treatment + 1 + (1|subject), data=clinic_info, method='PB')
summary(mod)
eta_squared(mod, partial = T)

library(emmeans)
mod <- afex::lmer_alt(HAMD ~ group + treatment + group*treatment + 1 + (1|subject), data=clinic_info, method='PB')
emm <- emmeans(mod, ~ group + treatment + group*treatment + 1)
contrast(emm, method = 'pairwise', adjust = 'fdr')
