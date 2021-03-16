library(ggRandomForests)
library("survival")
library("survminer")
suppressMessages(library(rpart))
suppressMessages(library(ggplot2))
suppressMessages(library(mlbench))
suppressMessages(library(caret))
suppressMessages(library(party))
suppressMessages(library(pROC))
suppressMessages(library(pROC))
library(randomForestSRC)
library(survival)
library(ranger)
library(ggplot2)
library(dplyr)
library(ggfortify)
library(ggRandomForests)

Data <- read.csv("/Users/jadonzhou/Research Projects/Healthcare Predictives 1/Variable variability studies/DM BMC Endocrine Disorders/Database.csv")
str(Data)

gg_dta <- gg_survival(interval = "Time", censor = "Event",data =Data,conf.int = 0.95)

# Kaplan Meier Survival Curve
km <- with(Data, Surv(Time, Event))
head(km,80)

dev.off()
plot(gg_dta) +labs(title='Kaplan-Meier survival curve: Overall cohort', y = "Survival Probability", x = "Age (Year)") +theme(legend.position = c(0.2, 0.2)) + coord_cartesian(y = c(0, 1.01))

ImData <- impute.rfsrc(data = Data)
# random survival forests
rfsrc <- rfsrc(Surv(Time, Event) ~ ., data = ImData,tree.err = TRUE,importance = TRUE)
plot(gg_error(rfsrc))

# variable importance
plot(gg_vimp(rfsrc)) + theme(legend.position = c(0.8, 0.2)) +labs(fill = "Variable importance > 0")+labs(title='Ophthalmological', y = "Variable importance", x = "Predictors")
plot(gg_vimp(rfsrc)) + theme(legend.position = c(0.8, 0.2)) +labs(fill = "Variable importance > 0")+labs(title='Ischemic stroke', y = "Variable importance", x = "Predictors")
plot(gg_vimp(rfsrc)) + theme(legend.position = c(0.8, 0.2)) +labs(fill = "Variable importance > 0")+labs(title='AF', y = "Variable importance", x = "Predictors")
plot(gg_vimp(rfsrc)) + theme(legend.position = c(0.8, 0.2)) +labs(fill = "Variable importance > 0")+labs(title='HF', y = "Variable importance", x = "Predictors")

varsel_pbc <- var.select(rfsrc)     
gg_md <- gg_minimal_depth(varsel_pbc)
print(gg_md)   
plot(gg_md)+labs(title='DM cohort', y = "Minimal depth", x = "Predictors")

plot(gg_minimal_vimp(gg_md)) + theme(legend.position=c(0.8, 0.2))+labs(title='(b) Overall cohort without interventions', y = "Variable importance ranking", x = "Minimal depth ranking order")


ggRFsrc <- plot(gg_rfsrc(rfsrc), alpha = 0.2) + theme(legend.position = "none") + labs(title='Predicted DM survivals with RSF model',y = "Survival Probability", x = "Age (Year)") + coord_cartesian(ylim = c(-0.01, 1.01))
show(ggRFsrc)

plot(gg_rfsrc(rfsrc, by = "Male.sex"))+
  theme(legend.position = c(0.2, 0.2))+
  labs(title='(b) Difference of predicted cardiac survivals between males and females',y = "Survival Probability", x = "Time",fill = "Sex")+
  coord_cartesian(ylim = c(-0.01, 1.01))


gg_v <- gg_variable(rfsrc, time = c(1000, 3000), time.labels = c("Observation time at 1000", "Observation time at 3000"))
plot(gg_v, xvar = "Age", alpha = 0.4) + theme(legend.position = "none") + coord_cartesian(ylim = c(-0.01, 1.01)) 


xvar.cat <- c("baseline_af")
plot(gg_v, xvar = xvar.cat, alpha = 0.4) + labs(y = "Survival") + theme(legend.position = "none") + coord_cartesian(ylim = c(-0.01, 1.02))

xvar.cat <- c("baseline_nyha", "baseline_af")
plot(gg_v, xvar = xvar.cat, alpha = 0.4) + labs(y = "Survival") + theme(legend.position = "none") + coord_cartesian(ylim = c(-0.01, 1.02))

# variable interactions
ggint <- gg_interaction(rfsrc)
plot(ggint, xvar = xvar)

partial_coplot_pbc <- gg_partial_coplot(rfsrc, xvar = "baseline_age", groups = ggvar$baseline_age, surv_type = "surv", time = rfsrc$time.interest[time_index[1000]], show.plots = FALSE)


R> ggplot(partial_coplot_pbc, aes(x=bili, y=yhat, col=group, shape=group)) +
  R+ geom_smooth(se = FALSE) +
  R+ labs(x = st.labs["bili"], y = "Survival at 1 year (%)",
          R+ color = "albumin", shape = "albumin")


ImData <- impute.rfsrc(data = Data, nodesize = 1, nsplit = 100)
r_fit <- ranger(Surv(Time, Event) ~ ., 
                data = ImData,
                importance = "permutation",
                splitrule = "extratrees",
                verbose = TRUE)
r_fit
# variable importance
vi <- data.frame(sort(round(r_fit$variable.importance, 4), decreasing = TRUE))
names(vi) <- "importance"
head(vi)

# Harrell's c-index
cat("Prediction Error = 1 - Harrell's c-index = ", r_fit$prediction.error)

v.obj1 <- rfsrc(Surv(Time, Event) ~ . , data = Data, importance = TRUE)
v.obj1
plot(v.obj1)

plot.survival(v.obj1, cens.model = "rfsrc")
plot.survival(v.obj1, subset = 3)

# Breiman-Cutler permutation vimp.
print(vimp(v.obj1)$importance)

# Breiman-Cutler random daughter vimp.
print(vimp(v.obj1, importance = "random")$importance)

# Breiman-Cutler joint permutation vimp.
print(vimp(v.obj1, joint = TRUE)$importance)


# visualization
v.obj <- rfsrc(Surv(Time, Event) ~ . , data = ImData, block.size = 1)
v.obj
plot(v.obj)




# For a decision tree, or random forest, scaling of data is not required.

##Run model

# Model will be run with Survived ~ Title.code + Group.code + Pclass + Sex + IsChildC12 + Embarked
#{r model-run}
# run model
set.seed(1234)

model <- cforest(Event ~ Age + Male + Mean.Fasting.Blood.Glucose + Mean.HbA1c + Baseline.Anemia + Total.protein,
                 data=train,
                 controls=cforest_unbiased(ntree=1500,mtry=4))

model

#

#{r fit-results}
# 
fitted <- predict(model, subset(train,select=c(3,4,5,6,7,8)), OOB=TRUE, type= "response")
levels(fitted) = c(0,1)

misClasificError <- mean(fitted != train$Survived)
print(paste('Training Accuracy', 1 - misClasificError))
#

# Model summary using #cforest# with 1500 trees
#{r print-model}

print(model)
#

# The ROC curve shows a value of 81.08%
#{r ROC-curve}
# ROC curve
predicted <- predict(model, train, OOB=TRUE, type= "response")
auc(as.numeric(train$Survived), as.numeric(predicted)  )
#

# The In-sample error shows 143 values miscatagorized with an accuracy of 83.95%.
#{r ein-error}
#
# compute in-sample results
caret::confusionMatrix(fitted, train$Survived)
#

# Plot of feature importance

#In order of importance Title.code, Sex and Pclass are the top 3 contributors. IsChildC12 contributes the least but it did add to the public score.

#{r importance-plot}
#plot feature importance
cforestImpPlot <- function(x) {
  cforest_importance <<- v <- varimp(x)
  dotchart(v[order(v)])
}

cforestImpPlot(model)
