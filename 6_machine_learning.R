library(data.table)
library(pROC)
library(ggplot2)

file = 'mechine_learning/total_before_after.csv'
df_org = fread(file)

compute_auc <- function(x, y){
  full_model = glm(y ~ x, family = binomial)
  probabilities <- predict(full_model, type = "response")
  roc_curve <- roc(y, probabilities)
  roc_data <- data.frame(
    FPR = roc_curve$specificities,  # False Positive Rate
    TPR = roc_curve$sensitivities,  # True Positive Rate
    thresholds = roc_curve$thresholds  # Thresholds
  )
  return(roc_data)
}
################################
df_before_org = df_org[treatment=='before', ]
df_before = df_before_org[,3:16]
X_before = scale(df_before[,2:14])
y_before = ifelse(df_before$group == 'patient', 1, 0)

################################
df_after_org = df_org[treatment=='after', ]
df_after = df_after_org[, c('group',"func_193_degree","func_192_dist","func_193_dist")]
X_after = scale(df_after[,2:4])
y_after = ifelse(df_after$group == 'patient', 1, 0)
################################
# refer to DBS group
# two groups: blue #8cb5c1, red f1afa0
################################
rocobj <- plot.roc(y_before, before_predict,
                   main="Before Tai Chi Exercise", percent=TRUE,
                   ci=TRUE, # compute AUC (of AUC by default)
                   print.auc=TRUE, asp = NA,
                   identity.col="darkgrey",
                   identity.lty=3,
                   identity.lwd=3,
) # print the AUC (will contain the CI)
ciobj <- ci.se(rocobj, # CI of sensitivity
               specificities=seq(0, 100, 5)) # over a select set of specificities
plot(ciobj, type="shape", col="#f5acc5") # plot as a blue shape
plot(ci.sp(rocobj, sensitivities=seq(0, 100, 5)), # ci of specificity
     xlim = c(0,100), ylim = c(0,100),
     type="bars") # print this one as bars




#### model comparison
model_before = glm(y_before ~ X_before, family = binomial)
model_after = glm(y_after ~ X_after, family = binomial)
before_predict = predict(model_before, type = "response")
after_predict = predict(model_after, type = "response")
rocobj1 <- plot.roc(y_before, before_predict,
                    main="Statistical comparison", percent=TRUE, col="#f5acc5", asp = NA)
rocobj2 <- lines.roc(y_after, after_predict, percent=TRUE, col="#f9f0bd")
#testobj <- roc.test(rocobj1, rocobj2)
#text(50, 50, labels=paste("p-value =", format.pval(testobj$p.value)), adj=c(0, .5))
#legend("bottomright", legend=c("Before", "After"), col=c("#1c61b6", "#008600"), lwd=2)



