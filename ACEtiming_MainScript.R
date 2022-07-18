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

source(here("functions", "randomForest_functions.R"))

######################## 
# load & inspect data
df <- read.csv2(here("data", "mainData_v1.csv"), dec = ".")
df$Group <- as.factor(df$Group)

str(df)
describe(df)
boxplot(select(df, starts_with("con")))



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

# all continuous data
corrM <- cor(df[, which(names(df) == "KERF_Sum"):ncol(df)], use = "pairwise.complete.obs")
corrplot(corrM)

# amygdala activity only
corrAmy <- cor(select(df, starts_with("con")))
corrplot.mixed(corrAmy)

# overall trauma, psychopathology and relevant amygdala activity
df_small1 <- df[, c("KERF_Sum", "CTQ_Sum", "BSI_GSI", "PCL_5_Sum", "FDS_Sum", 
                    "con0009_AAL_Amygdala_right", "con0009_AAL_Amygdala_left", "con0010_AAL_Amygdala_right", "con0010_AAL_Amygdala_left")]
corrSelect <- cor(df_small1)
correSelect_pAdjust <- corr.p(corrSelect, n = nrow(df_small1), adjust = "holm")$p
corrplot(corrSelect, p.mat = correSelect_pAdjust) # below diagonal: uncorrected p-values. above diagonal: corrected p-values 

plot(df$KERF_Sum, df$con0010_AAL_Amygdala_right)
cor.test(df$KERF_Sum, df$con0010_AAL_Amygdala_right)

plot(df$FDS_Sum, df$con0009_AAL_Amygdala_right)
cor.test(df$FDS_Sum, df$con0009_AAL_Amygdala_right)

plot(df$FDS_Sum, df$con0009_AAL_Amygdala_left)
cor.test(df$FDS_Sum, df$con0009_AAL_Amygdala_left)



########################
# random forests
ACErandomForest(DV = "con0010_AAL_Amygdala_right", include.imp = F)
ACErandomForest(DV = "con0010_AAL_Amygdala_left", include.imp = F)
ACErandomForest(DV = "con0009_AAL_Amygdala_right", include.imp = F)
ACErandomForest(DV = "con0009_AAL_Amygdala_left", include.imp = F)

########################
# permutation loops
DVs <- c("con0010_AAL_Amygdala_right", "con0010_AAL_Amygdala_left", "con0009_AAL_Amygdala_right", "con0009_AAL_Amygdala_left")

for(i in 1:length(DVs)){
  RF.permutation(DV = DVs[i])
}

