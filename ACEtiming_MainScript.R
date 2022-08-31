########################
# TO DO

# repeated random forest for better accuracy (maybe only for main estimate)
# run permutation test
# check out relevance of mincriterion
# model comparison when time etc. is included



######################## 
# load packages & functions

library("here")
library("psych")
library("corrplot")
library("tidyselect")
library("dplyr")
library("reshape2")
library("ggplot2")
library("party")
library("Metrics")
library("psychometric")
library("tidyverse")
library("gdata")
library("mediation")
library("ez")

source(here("functions", "randomForest_functions.R"))

######################## 
# load & inspect data
df <- read.csv2(here("data", "mainData_v2.csv"), dec = ".")
df$Group <- as.factor(df$Group)

options(max.print = 4000)

str(df)
describe(df)
boxplot(dplyr::select(df, starts_with("AAL")))


########################
# check amygdala reactivity to different emotions

apply(dplyr::select(df, starts_with("AAL")), 2, t.test)

dfANOVA <- df %>% 
  dplyr::select("ID", starts_with("AAL")) %>% 
  pivot_longer(starts_with("AAL"))

dfANOVA <- data.frame(dfANOVA, str_split_fixed(dfANOVA$name, "_", 6))
  
amyModel <- ezANOVA(data = dfANOVA, 
                    dv = value,
                    within = .(X3, X4),
                    wid = ID)
print(amyModel)

######################## 
# plot ACE

df_timing <- data.frame(df$ID, df[, which(names(df) == "KERF_Sum_0"):which(names(df) == "KERF_Sum_17")])
df_timing_long <- melt(data = df_timing, id.vars = "df.ID")

ggplot(data = df_timing_long, aes(x = variable, y = value, group = df.ID)) +
  #draw individual data
  geom_line(alpha=0.12, size = 1) +
  # create white background for mean data
  stat_summary(aes(group = 1), geom = "line", fun = mean, size = 1.2, colour = "white") + 
  stat_summary(aes(group = 1), geom = "point", fun = mean, size = 4.2, colour = "white") + 
  stat_summary(aes(group = 1), fun.data = mean_se, fun.args = list(mult = 1.96), geom = "errorbar", width = 0.55, size = 1.2, colour = "white") +
  # draw mean values on white background
  stat_summary(aes(group = 1), geom = "line", fun = mean, size = 1) +
  stat_summary(aes(group = 1), geom = "point", fun = mean, size = 4) +
  stat_summary(aes(group = 1), fun.data = mean_se, fun.args = list(mult = 1.96), geom = "errorbar", width = 0.5, size = 1) +

  xlab("Age at adversity") + ylab("MACE") + 
  labs(color = NULL) +
  scale_colour_discrete(labels = c("PTSD", "Trauma Controls"), guide = guide_legend(reverse=FALSE)) +
  scale_x_discrete(labels= as.factor(seq(from = 0, to = 17, by = 1))) +
  
  theme_classic() +
  theme(legend.position = "top")



########################
# correlations

# amygdala activity only
corrAmy <- cor(dplyr::select(df, starts_with("AAL")))
corrplot.mixed(corrAmy)

# overall trauma, psychopathology and relevant amygdala activity
df_small1 <- dplyr::select(df, c("KERF_Sum", "CTQ_Sum", "BSI_GSI", "PCL_5_Sum", "FDS_ges_Sum", starts_with("AAL")))

corrSelect <- cor(df_small1)
correSelect_pAdjust <- corr.p(corrSelect, n = nrow(df_small1), adjust = "holm")$p
corrplot(corrSelect, p.mat = correSelect_pAdjust) # below diagonal: uncorrected p-values. above diagonal: corrected p-values 

plot(df$KERF_Sum, df$con0009_AAL_Amygdala_right_anger_fear_vs_shapes)
cor.test(df$KERF_Sum, df$con0009_AAL_Amygdala_right_anger_fear_vs_shapes)

plot(df$BSI_GSI, df$con0009_AAL_Amygdala_right_anger_fear_vs_shapes)
cor.test(df$BSI_GSI, df$con0009_AAL_Amygdala_right_anger_fear_vs_shapes)

plot(df$FDS_ges_Sum, df$con0009_AAL_Amygdala_right_anger_fear_vs_shapes)
cor.test(df$FDS_ges_Sum, df$con0009_AAL_Amygdala_right_anger_fear_vs_shapes)



########################
# random forests for anger&fear vs shapes
ACErandomForest(DV = "con0009_AAL_Amygdala_right_anger_fear_vs_shapes", include.imp = F)
ACErandomForest(DV = "con0009_AAL_Amygdala_left_anger_fear_vs_shapes", include.imp = F)

# reapeat and average results for right amygdala for stability
varExpl <- numeric(10)
for(i in 1:10){
  set.seed(777+i)
  varExpl[i] <- ACErandomForest(DV = "con0009_AAL_Amygdala_right_anger_fear_vs_shapes", include.imp = F)$Variance_explained
}
mean(varExpl)

# check for contrast anger versus shapes
varExpl <- numeric(10)
for(i in 1:10){
  set.seed(888+i)
  varExpl[i] <- ACErandomForest(DV = "AAL_Amygdala_right_angry_vs_shapes", include.imp = F)$Variance_explained
}
mean(varExpl)

# check for contrast fearful versus shapes
varExpl <- numeric(10)
for(i in 1:10){
  set.seed(888+i)
  varExpl[i] <- ACErandomForest(DV = "AAL_Amygdala_right_fearful_vs_shapes", include.imp = F)$Variance_explained
}
mean(varExpl)


ACErandomForest(DV = "BSI_GSI", include.imp = F)
ACErandomForest(DV = "FDS_ges_Sum", include.imp = F)



########################
# model comparisons
#[at the moment with all life years for the typeTiming model!]

model1_baseline <- c("KERF_Sum", "KERF_Multi", "KERF_Duration", "Group", "Age", "Sex")
model2_type <- c(model1_baseline, names(df)[which(names(df) == "KERF_SUM_PEA"):which(names(df) == "KERF_SUM_SEXA_O")])
model3_timing <- c(model2_type, names(df)[which(names(df) == "KERF_Sum_3"):which(names(df) == "KERF_Sum_17")])
model4_typeTiming <- c(model2_type, names(df)[which(names(df) == "SUM_PEA_0"):which(names(df) == "SUM_WITS_17")])
model5_psychoPath <- c(model1_baseline, "BSI_GSI", "BDI.II_Sum", "PCL_5_Sum", "SSD.12_Sum", "FDS_ges_Sum", "FDS_DES_Sum")
model6_psychoPatTypeTiming <- c(model5_psychoPath, names(df)[which(names(df) == "SUM_PEA_0"):which(names(df) == "SUM_WITS_17")])

modelList <- list(model1_baseline, model2_type, model3_timing, model4_typeTiming, model5_psychoPath, model6_psychoPatTypeTiming)


Rf.compareModels("con0009_AAL_Amygdala_right_anger_fear_vs_shapes", 
                 predictorSets = modelList,
                 mtryArg = "sqroot")

Rf.compareModels("con0009_AAL_Amygdala_left_anger_fear_vs_shapes", 
                 predictorSets = modelList,
                 mtryArg = "sqroot")

Rf.compareModels("BSI_GSI", 
                 predictorSets = modelList,
                 mtryArg = "sqroot")

psychoPathModel <- cforest(data = df[, c("con0009_AAL_Amygdala_right_anger_fear_vs_shapes", model5_psychoPath)],
                              con0009_AAL_Amygdala_right_anger_fear_vs_shapes ~ .,
                              controls = cforest_unbiased(mtry = sqrt(length(model5_psychoPath))))

plot(varimp(psychoPathModel, conditional = FALSE))


########################
# mediation

fit.totaleffect <- lm(data = df,
                   con0009_AAL_Amygdala_right_anger_fear_vs_shapes ~ KERF_Sum)
summary(fit.totaleffect)

fit.mediator <- lm(data = df,
                FDS_ges_Sum ~ KERF_Sum)
summary(fit.mediator)

fit.dv <- lm(data = df,
             con0009_AAL_Amygdala_right_anger_fear_vs_shapes ~ FDS_ges_Sum+KERF_Sum)
summary(fit.dv)

medResults = mediation::mediate(fit.mediator, fit.dv, treat='KERF_Sum', mediator='FDS_ges_Sum', boot=T)
summary(medResults)



interactionModel <- lm(data = df,
                        con0009_AAL_Amygdala_right_anger_fear_vs_shapes ~ FDS_ges_Sum*KERF_Sum)
summary(interactionModel)

