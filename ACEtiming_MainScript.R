########################
# TO DO

# repeated random forest for better accuracy (maybe only for main estimate)
# run permutation test
# correlation replications of previous studies
# remove group variable from all models 
# adapt rf.compare function so that accuracy vector is named with model names if multiple models are provided


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
library("Hmisc")
library("gameofthrones")
library("cocor")
library("cowplot")

source(here("functions", "randomForest_functions.R"))

######################## 
# load & inspect data
df <- read.csv2(here("data", "mainData_v3.csv"), dec = ".")
df$Group <- as.factor(df$Group)

options(max.print = 4000)

str(df)
describe(df)
boxplot(dplyr::select(df, starts_with("AAL")))



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
  
  xlab("Age at Childhood Maltreatment Exposure") + ylab("KERF-40+ Sum Score") + 
  labs(color = NULL) +
  scale_colour_discrete(labels = c("PTSD", "Trauma Controls"), guide = guide_legend(reverse=FALSE)) +
  scale_x_discrete(labels= as.factor(seq(from = 0, to = 17, by = 1))) +
  
  theme_classic() +
  theme(legend.position = "top")

ggsave(here("figures", "Figure1_ACEexposureTimings.pdf"), device = "pdf")
ggsave(here("figures", "Figure1_ACEexposureTimings.png"), device = "tiff")



########################
# check amygdala reactivity to different emotions

df_amyOnly <- df[, c("con0009_AAL_Amygdala_right_anger_fear_vs_shapes", "con0009_AAL_Amygdala_left_anger_fear_vs_shapes", 
                     "con0012_AAL_Amygdala_right_neutral_surprised_vs_shapes", "con0012_AAL_Amygdala_left_neutral_surprised_vs_shapes")]

apply(df_amyOnly, 2, t.test)
apply(df_amyOnly, 2, function(x) mean(x)/sd(x))

dfANOVA <- melt(data.frame(ID = df$ID, df_amyOnly), id.vars = "ID")
dfANOVA <- data.frame(dfANOVA, str_split_fixed(dfANOVA$variable, "_", 6))
names(dfANOVA)[c(7,9)] <- c("hemisphere", "valence")
  
amyModel <- ezANOVA(data = dfANOVA, 
                    dv = value,
                    within = .(hemisphere, valence),
                    wid = ID)
print(amyModel)

rcorr(as.matrix(df_amyOnly))


########################
# model comparisons
# sources: Yeh, A. (2000). More accurate tests for the statistical significance of result differences. arXiv preprint cs/0008005
# https://cborchers.com/2021/04/22/permutation-test-for-f-score-differences-in-python/

model1_baseline <- c("KERF_Sum", "KERF_Multi", "KERF_Duration", "Age", "Sex", "Scanner", "Medication_Load")
model2_type <- c(model1_baseline, names(df)[which(names(df) == "KERF_SUM_PEA"):which(names(df) == "KERF_SUM_SEXA_O")])
model3_timing <- c(model1_baseline, names(df)[which(names(df) == "KERF_Sum_2"):which(names(df) == "KERF_Sum_17")])
model4_typeTiming <- c(model2_type, names(df)[which(names(df) == "SUM_PEA_0"):which(names(df) == "SUM_WITS_17")]) 
model4_typeTiming <- model4_typeTiming[str_detect(model4_typeTiming, "_0|_1$", negate = TRUE)] # 2nd line for model 4 removes ages 0-1
model5_psychoPath <- c(model1_baseline, "BSI_GSI", "BDI.II_Sum", "PCL_5_Sum", "SSD.12_Sum", "FDS_DES_Sum")
model6_psychoPatTypeTiming <- unique(c(model5_psychoPath, model4_typeTiming))

modelList <- list(model1_baseline, model2_type, model3_timing, model4_typeTiming, model5_psychoPath, model6_psychoPatTypeTiming)

# make dataframe without missings on any predictors of interest for any model of interest
dfcomplete <- df[complete.cases(df[,model6_psychoPatTypeTiming]), ]


############
# random forest regression: Descriptive variance explained


# predict right amygdala [might take a couple of minutes]
amyRight_AccRepeat <- RFmain_repeat_multiModel("con0009_AAL_Amygdala_right_anger_fear_vs_shapes", 
                 predictorSets = modelList,
                 mtryArg = "bagging",
                 data = dfcomplete,
                 repeats = 5,
                 seed = 1000)

# predict left amygdala
amyLeft_AccRepeat <- RFmain_repeat_multiModel("con0009_AAL_Amygdala_left_anger_fear_vs_shapes", 
                                               predictorSets = modelList,
                                               mtryArg = "bagging",
                                               data = dfcomplete,
                                               repeats = 5,
                                               seed = 1001)

# predict right amygdala of control contrast
amyRight_control_AccRepeat <- RFmain_repeat_multiModel("con0012_AAL_Amygdala_right_neutral_surprised_vs_shapes", 
                                               predictorSets = modelList,
                                               mtryArg = "bagging",
                                               data = dfcomplete,
                                               repeats = 5,
                                               seed = 1002)

# predict left amygdala of control contrast
amyLeft_control_AccRepeat <- RFmain_repeat_multiModel("con0012_AAL_Amygdala_left_neutral_surprised_vs_shapes", 
                                              predictorSets = modelList,
                                              mtryArg = "bagging",
                                              data = dfcomplete,
                                              repeats = 5,
                                              seed = 1003)

# check results
amyRight_AccRepeat
amyLeft_AccRepeat
amyRight_control_AccRepeat
amyLeft_control_AccRepeat


############
# create permutations for significance tests
# [commented out, because it takes ~10-40h per model to do on our hardware]

# amyRight_model5_perm <- RFperm("con0009_AAL_Amygdala_right_anger_fear_vs_shapes",
#                                predictorSets = modelList[[5]],
#                                data = dfcomplete,
#                                mtryArg = "bagging",
#                                numtree = 1000,
#                                nPerm = 1000)
# 
# if(!file.exists(here("results", "amyRight_model5_perm.csv"))){
#   write.csv(amyRight_model5_perm, here("results", "amyRight_model5_perm.csv"), row.names = FALSE)
# }else{
#   warning("file already exists in folder")
# }


# amyLeft_model5_perm <- RFperm("con0009_AAL_Amygdala_left_anger_fear_vs_shapes",
#                                predictorSets = modelList[[5]],
#                                data = dfcomplete,
#                                mtryArg = "bagging",
#                                numtree = 1000,
#                                nPerm = 1000)
# 
# if(!file.exists(here("results", "amyLeft_model5_perm.csv"))){
#   write.csv(amyLeft_model5_perm, here("results", "amyLeft_model5_perm.csv"), row.names = FALSE)
# }else{
#   warning("file already exists in folder")
# }


# amyLeft_model6_perm <- RFperm("con0009_AAL_Amygdala_left_anger_fear_vs_shapes",
#                                predictorSets = modelList[[6]],
#                                data = dfcomplete,
#                                mtryArg = "bagging",
#                                numtree = 1000,
#                                nPerm = 1000)
# 
# if(!file.exists(here("results", "amyLeft_model6_perm.csv"))){
#   write.csv(amyLeft_model6_perm, here("results", "amyLeft_model6_perm.csv"), row.names = FALSE)
# }else{
#   warning("file already exists in folder")
# }



############
# calculate p-values of model 5 for both hemispheres

# right amygdala
RFrep_right_model5 <- RFmain_repeat_singleModel("con0009_AAL_Amygdala_right_anger_fear_vs_shapes", 
                                                predictorSets = modelList[[5]],
                                                data = dfcomplete,
                                                mtryArg = "bagging",
                                                include.imp = TRUE,
                                                repeats = 5,
                                                seed = 666)

model5Perms_right <- read.csv(here("results", "amyRight_model5_perm.csv"))
cbind(names(model5Perms_right), returnP(c(amyRight_AccRepeat[5], RFrep_right_model5$Importance), model5Perms_right))
amyRight_model5_pValues <- cbind(names(model5Perms_right), returnP(c(amyRight_AccRepeat[5], RFrep_right_model5$Importance), model5Perms_right))

if(!file.exists(here("results", "amyRight_model5_pValues.csv"))){
  write.csv(amyRight_model5_pValues, here("results", "amyRight_model5_pValues.csv"), row.names = FALSE)
}else{
  warning("file already exists in folder")
}


# left amygdala
RFrep_left_model5 <- RFmain_repeat_singleModel("con0009_AAL_Amygdala_left_anger_fear_vs_shapes", 
                                                predictorSets = modelList[[5]],
                                                data = dfcomplete,
                                                mtryArg = "bagging",
                                                include.imp = TRUE,
                                                repeats = 5,
                                                seed = 666)

model5Perms_left <- read.csv(here("results", "amyLeft_model5_perm.csv"))
cbind(names(model5Perms_left), returnP(c(amyLeft_AccRepeat[5], RFrep_left_model5$Importance), model5Perms_left))
amyLeft_model5_pValues <- cbind(names(model5Perms_left), returnP(c(amyLeft_AccRepeat[5], RFrep_left_model5$Importance), model5Perms_left))

if(!file.exists(here("results", "amyLeft_model5_pValues.csv"))){
  write.csv(amyLeft_model5_pValues, here("results", "amyLeft_model5_pValues.csv"), row.names = FALSE)
}else{
  warning("file already exists in folder")
}


############
# calculate p-values of model 6 for both hemispheres


# right amygdala
RFrep_right_model6 <- RFmain_repeat_singleModel("con0009_AAL_Amygdala_right_anger_fear_vs_shapes", 
                                                predictorSets = modelList[[6]],
                                                data = dfcomplete,
                                                mtryArg = "bagging",
                                                include.imp = TRUE,
                                                repeats = 5,
                                                seed = 666)

model6Perms_right <- read.csv(here("results", "amyRight_model6_perm.csv"))
cbind(names(model6Perms_right), returnP(c(amyRight_AccRepeat[6], RFrep_right_model6$Importance), model6Perms_right))
amyRight_model6_pValues <- cbind(names(model6Perms_right), returnP(c(amyRight_AccRepeat[6], RFrep_right_model6$Importance), model6Perms_right))

if(!file.exists(here("results", "amyRight_model6_pValues.csv"))){
  write.csv(amyRight_model6_pValues, here("results", "amyRight_model6_pValues.csv"), row.names = FALSE)
}else{
  warning("file already exists in folder")
}


# left amygdala
RFrep_left_model6 <- RFmain_repeat_singleModel("con0009_AAL_Amygdala_left_anger_fear_vs_shapes", 
                                               predictorSets = modelList[[6]],
                                               data = dfcomplete,
                                               mtryArg = "bagging",
                                               include.imp = TRUE,
                                               repeats = 5,
                                               seed = 666)

model6Perms_left <- read.csv(here("results", "amyLeft_model6_perm.csv"))
cbind(names(model6Perms_left), returnP(c(amyLeft_AccRepeat[6], RFrep_left_model6$Importance), model6Perms_left))
amyLeft_model6_pValues <- cbind(names(model6Perms_left), returnP(c(amyLeft_AccRepeat[6], RFrep_left_model6$Importance), model6Perms_left))

if(!file.exists(here("results", "amyLeft_model6_pValues.csv"))){
  write.csv(amyLeft_model6_pValues, here("results", "amyLeft_model6_pValues.csv"), row.names = FALSE)
}else{
  warning("file already exists in folder")
}



############
# model comparisons for significant models


# compare Variance explained of Models (set "saveModels" to TRUE to later compare models for significance)
amyRightModels <- RFmain("con0009_AAL_Amygdala_right_anger_fear_vs_shapes", 
                         predictorSets = modelList,
                         mtryArg = "bagging",
                         data = dfcomplete,
                         saveModels = TRUE)
amyRightModels$Accuracies

# compare two models for differences in accuracy
# [using the accuracies from the more stable repeated models]
RFcompare(Diff = (amyRight_AccRepeat[5] - amyRight_AccRepeat[6]),
          responseVar = dfcomplete$con0009_AAL_Amygdala_right_anger_fear_vs_shapes,
          model1 = amyRightModels$RFmodels[[5]],
          model2 = amyRightModels$RFmodels[[6]],
          nPerms = 1000)

# left amygdala model comparison
amyLeftModels <- RFmain("con0009_AAL_Amygdala_left_anger_fear_vs_shapes", 
                        predictorSets = modelList,
                        mtryArg = "bagging",
                        data = dfcomplete,
                        saveModels = TRUE)
amyLeftModels$Accuracies

RFcompare(Diff = (amyLeft_AccRepeat[5] - amyLeft_AccRepeat[6]),
          responseVar = dfcomplete$con0009_AAL_Amygdala_left_anger_fear_vs_shapes,
          model1 = amyLeftModels$RFmodels[[5]],
          model2 = amyLeftModels$RFmodels[[6]],
          nPerms = 1000)




############
# plot variance explained for different models

modelNames <- c("Model 1: Baseline", "Model 2: Type", "Model 3: Timing", "Model 4: Type x Timing", "Model 5: Psychopathology", "Model 6: Full Model")
df_plotAcc <- data.frame("modelNames" = rep(modelNames, 2), "Accuracies" = c(amyRight_AccRepeat, amyLeft_AccRepeat), "Hemisphere" = rep(c("Right", "Left"), each = length(modelNames)))

ggplot(data = df_plotAcc, aes(y = Accuracies, x = modelNames, fill = Hemisphere)) +
  geom_bar(position = "dodge", stat = "identity", colour = "black") +
  
  ylim(-0.2, 0.4) + 
  ylab(expression(paste(italic("R²"), "(cross-validated)"))) +

  xlab(NULL) +
  
  theme_classic() +
  scale_fill_got_d(option = "Daenerys") +
  theme(axis.text.x=element_text(angle = 45, vjust = 0.5)) +
  
  annotate("text", x = 4.775, y = 0.035, label = "*") +
  annotate("text", x = 5.225, y = 0.075, label = "***") +
  annotate("text", x = 5.775, y = 0.06, label = "**") +
  annotate("text", x = 6.225, y = 0.10, label = "**") +
  
  geom_segment(aes(x = 4.66, y = 0.15, xend = 5.85, yend = 0.15)) +
  annotate("text", x = 5.25, y = 0.165, label = "N.S.") +
  
  geom_segment(aes(x = 5.1, y = 0.22, xend = 6.3, yend = 0.22)) +
  annotate("text", x = 5.75, y = 0.235, label = "N.S.") +
  
  scale_y_continuous(n.breaks = 12, limits = c(-0.2, 0.4))
  
ggsave(here("figures", "Figure2_AccuracyModelComparisons.png"), device = "png")
ggsave(here("figures", "Figure2_AccuracyModelComparisons.pdf"), device = "pdf")
  



############
# plot variable importance for model 5



predictorList <- names(RFrep_right_model5$Importance)
df_plotImp <- data.frame("Predictors" = rep(predictorList, 2), "varImp" = c(RFrep_right_model5$Importance, RFrep_left_model5$Importance), "Hemisphere" = rep(c("Right", "Left"), each = length(predictorList)))
df_plotImp$Predictors <-  factor(df_plotImp$Predictors, unique(df_plotImp$Predictors))
df_plotImp$PredLabels <- c("KERF-40+ Sum", "KERF-40+ Multiplicity", "KERF-40+ Duration", "Age", "Sex", "MRI Scanner Type", "Medication Load", "BSI GSI", "BDI-II", "PCL-5", "SSD-12", "FDS DES")
df_plotImp$PredLabels <-  factor(df_plotImp$PredLabels, unique(df_plotImp$PredLabels))

ggplot(data = df_plotImp, aes(y = varImp, x = PredLabels, fill = Hemisphere)) +
  geom_bar(position = "dodge", stat = "identity", colour = "black") +
  
  ylab("Variable Importance") +
  
  xlab(NULL) +
  
  theme_classic() +
  scale_fill_got_d(option = "Daenerys") +
  theme(axis.text.x=element_text(angle = 45, vjust = 0.5)) +
  
  scale_y_continuous(n.breaks = 12, limits = c(-0.04, 0.2)) +
  annotate("text", x = 11.8, y = 0.14, label = "***") +
  annotate("text", x = 12.3, y = 0.115, label = "***")


ggsave(here("figures", "Figure3_variableImportance.png"), device = "png")
ggsave(here("figures", "Figure3_variableImportance.pdf"), device = "pdf")



########################
# bivariate association with FDS
cor.test(df$FDS_DES_Sum, df$con0009_AAL_Amygdala_right_anger_fear_vs_shapes)
cor.test(df$FDS_DES_Sum, df$con0009_AAL_Amygdala_left_anger_fear_vs_shapes)

lAmyPlot <- ggplot(data = df, aes(x = FDS_DES_Sum, y = con0009_AAL_Amygdala_left_anger_fear_vs_shapes)) + 
  geom_point(alpha = 0.3) + 
  geom_smooth(method = "lm", colour="black") +
  
  xlab(NULL) + 
  ylab("Contrast Values \n(Threatening Faces versus Shapes)") + 
  
  theme_classic() +
  
  ggtitle("Left Amygdala") +
  theme(plot.title = element_text(hjust = 0.5))
lAmyPlot

rAmyPlot <- ggplot(data = df, aes(x = FDS_DES_Sum, y = con0009_AAL_Amygdala_right_anger_fear_vs_shapes)) + 
  geom_point(alpha = 0.3) + 
  geom_smooth(method = "lm", colour="black") +
  
  xlab(NULL) + 
  #ylab("contrast values", colour = "white") + 
  
  theme_classic() +
  
  ggtitle("Right Amygdala") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.y = element_text(colour = "white"))

twoPanelFig <- plot_grid(lAmyPlot, rAmyPlot)
ggdraw(add_sub(twoPanelFig, "Trait Dissociation\n(FDS DES)", size = 12))

ggsave(here("figures", "Figure4_FDScorrelations.png"), device = "png")
ggsave(here("figures", "Figure4_FDScorrelations.pdf"), device = "pdf")




########################
# mediation

cor.test(df$FDS_DES_Sum, df$KERF_Sum)
cor.test(df$KERF_Sum, df$con0009_AAL_Amygdala_right_anger_fear_vs_shapes)
cor.test(df$KERF_Sum, df$con0009_AAL_Amygdala_left_anger_fear_vs_shapes)

# right amygdala
dfMediation <- df[, c("con0009_AAL_Amygdala_right_anger_fear_vs_shapes", "con0009_AAL_Amygdala_left_anger_fear_vs_shapes", "KERF_Sum", "KERF_Multi", "KERF_Duration", "FDS_DES_Sum")]
dfMediation <- data.frame(scale(dfMediation))

fit.totaleffect_right <- lm(data = dfMediation,
                   con0009_AAL_Amygdala_right_anger_fear_vs_shapes ~ KERF_Sum+Sex+Age)
summary(fit.totaleffect_right)

fit.mediator_right <- lm(data = dfMediation,
                FDS_DES_Sum ~ KERF_Sum+Sex+Age)
summary(fit.mediator_right)

fit.dv_right <- lm(data = dfMediation,
             con0009_AAL_Amygdala_right_anger_fear_vs_shapes ~ FDS_DES_Sum+KERF_Sum+Sex+Age)
summary(fit.dv_right)

medResults_right = mediation::mediate(fit.mediator_right, fit.dv_right, treat='KERF_Sum', mediator='FDS_DES_Sum', boot=T, sims = 5000)
summary(medResults_right)


# left amygdala
fit.totaleffect_left <- lm(data = dfMediation,
                            con0009_AAL_Amygdala_left_anger_fear_vs_shapes ~ KERF_Sum+Sex+Age)
summary(fit.totaleffect_left)

fit.mediator_left <- lm(data = dfMediation,
                         FDS_DES_Sum ~ KERF_Sum+Sex+Age)
summary(fit.mediator_left)

fit.dv_left <- lm(data = dfMediation,
                   con0009_AAL_Amygdala_left_anger_fear_vs_shapes ~ FDS_DES_Sum+KERF_Sum+Sex+Age)
summary(fit.dv_left)

medResults_left = mediation::mediate(fit.mediator_left, fit.dv_left, treat='KERF_Sum', mediator='FDS_DES_Sum', boot=T, sims = 5000)
summary(medResults_left)


########################
# random forests predicting dissociation

RF_FDS_models <- RFmain("FDS_DES_Sum",
       predictorSets = modelList[1:4],
       mtryArg = "bagging",
       data = dfcomplete,
       saveModels = TRUE)


RF_FDS_repAccs <- RFmain_repeat_multiModel("FDS_DES_Sum", 
                         predictorSets = modelList[1:4],
                         mtryArg = "bagging",
                         data = dfcomplete,
                         repeats = 5,
                         seed = 1001)

RF_FDS_repSingle <- RFmain_repeat_singleModel("FDS_DES_Sum", 
                                           predictorSets = modelList[[1]],
                                           mtryArg = "bagging",
                                           data = dfcomplete,
                                           repeats = 5,
                                           include.imp = TRUE,
                                           seed = 1001)


# create dataframe with model accuracy and variable importances from permuted model
# [save the results to another folder, because this procedure takes long]
FDS_model1_perm <- RFperm("FDS_DES_Sum",
                               predictorSets = modelList[[1]],
                               data = dfcomplete,
                               mtryArg = "bagging",
                               numtree = 1000,
                               nPerm = 1000)

if(!file.exists(here("results", "FDS_model1_perm.csv"))){
  write.csv(FDS_model1_perm, here("results", "FDS_model1_perm.csv"), row.names = FALSE)
}else{
  warning("file already exists in folder")
}

# compute p-values
FDS_model1_perm <- read.csv(here("results", "FDS_model1_perm.csv"))
cbind(names(FDS_model1_perm), returnP(c(RF_FDS_repAccs[1], RF_FDS_repSingle$Importance), FDS_model1_perm))

set.seed(100)
RFcompare(Diff = (RF_FDS_repAccs[1] - RF_FDS_repAccs[2]),
          responseVar = dfcomplete$FDS_DES_Sum,
          model1 = RF_FDS_models$RFmodels[[1]],
          model2 = RF_FDS_models$RFmodels[[2]],
          nPerms = 1000)

RF_FDS_repSingle_type <- RFmain_repeat_singleModel("FDS_DES_Sum", 
                                              predictorSets = modelList[[2]],
                                              mtryArg = "bagging",
                                              data = dfcomplete,
                                              repeats = 5,
                                              include.imp = TRUE,
                                              seed = 1001)
RF_FDS_repSingle_type


########################
# life years of interest

df$amyBilat_negFaces <- (df$con0009_AAL_Amygdala_right_anger_fear_vs_shapes + df$con0009_AAL_Amygdala_left_anger_fear_vs_shapes)/2

df$teicher_prePub_parentPeer <- (df$SUM_PPA_2+df$SUM_PPA_3+df$SUM_PEER_5)/3
df$teicher_postPub <- (df$SUM_PEER_12+df$SUM_PEER_14)/2

cor.test(df$teicher_prePub_parentPeer, df$amyBilat_negFaces)

cor.test(df$SUM_PPA_2, df$con0009_AAL_Amygdala_right_anger_fear_vs_shapes)
cor.test(df$SUM_PPA_2, df$con0009_AAL_Amygdala_left_anger_fear_vs_shapes)
cor.test(df$SUM_PPA_3, df$con0009_AAL_Amygdala_right_anger_fear_vs_shapes)
cor.test(df$SUM_PPA_3, df$con0009_AAL_Amygdala_left_anger_fear_vs_shapes)
cor.test(df$SUM_PPA_5, df$con0009_AAL_Amygdala_right_anger_fear_vs_shapes)
cor.test(df$SUM_PPA_5, df$con0009_AAL_Amygdala_left_anger_fear_vs_shapes)

cor.test(df$SUM_PEER_2, df$con0009_AAL_Amygdala_right_anger_fear_vs_shapes)
cor.test(df$SUM_PEER_2, df$con0009_AAL_Amygdala_left_anger_fear_vs_shapes)
cor.test(df$SUM_PEER_3, df$con0009_AAL_Amygdala_right_anger_fear_vs_shapes)
cor.test(df$SUM_PEER_3, df$con0009_AAL_Amygdala_left_anger_fear_vs_shapes)
cor.test(df$SUM_PEER_5, df$con0009_AAL_Amygdala_right_anger_fear_vs_shapes)
cor.test(df$SUM_PEER_5, df$con0009_AAL_Amygdala_left_anger_fear_vs_shapes)

cor.test(df$teicher_postPub, df$amyBilat_negFaces)

cor.test(df$SUM_PEER_12, df$con0009_AAL_Amygdala_right_anger_fear_vs_shapes)
cor.test(df$SUM_PEER_12, df$con0009_AAL_Amygdala_left_anger_fear_vs_shapes)
cor.test(df$SUM_PEER_14, df$con0009_AAL_Amygdala_right_anger_fear_vs_shapes)
cor.test(df$SUM_PEER_14, df$con0009_AAL_Amygdala_left_anger_fear_vs_shapes)


df$sicorello_prePub <- (df$KERF_Sum_3+df$KERF_Sum_4)/2
df$sicorello_postPub <- (df$KERF_Sum_16+df$KERF_Sum_17)/2

summary(
  lm(data = df,
     scale(con0009_AAL_Amygdala_right_anger_fear_vs_shapes) ~ scale(sicorello_prePub)+scale(PCL_5_Sum)
  )
)

summary(
  lm(data = df,
     scale(con0009_AAL_Amygdala_right_anger_fear_vs_shapes) ~ scale(sicorello_postPub)+scale(PCL_5_Sum)
  )
)

summary(
  lm(data = df,
     scale(con0009_AAL_Amygdala_right_anger_fear_vs_shapes) ~ scale(sicorello_prePub)*scale(PCL_5_Sum)
  )
)

summary(
  lm(data = df,
     scale(con0009_AAL_Amygdala_right_anger_fear_vs_shapes) ~ scale(sicorello_postPub)*scale(PCL_5_Sum)
  )
)

cor.test(df$sicorello_prePub, df$con0009_AAL_Amygdala_right_anger_fear_vs_shapes)
cor.test(df$sicorello_postPub, df$con0009_AAL_Amygdala_right_anger_fear_vs_shapes)
