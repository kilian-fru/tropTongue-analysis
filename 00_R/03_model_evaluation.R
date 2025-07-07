library(tidyverse)
library(BeeIT)
library(MuMIn)
library(gridExtra)
library(grid)
library(lme4)

rm(list=ls())

##### import trait data & create dataframe for model training #####

traits <- read.csv2("02_processed_data/traits.csv")

prediction <- traits %>%
  pivot_wider(id_cols=individual:species, names_from=trait, 
              values_from=value)

##### create linear model for Meliponini #####

# pred_mel_training <- prediction %>% 
#   filter(tribe=="Meliponini") %>%
#   sample_n(0.8*343) %>%
#   group_by(tribe, species) %>%
#   mutate(numb = n()) %>%
#   ungroup()
# 
# write.csv2(pred_mel_training, "02_processed_data/training_test_model.csv")

pred_mel_training <- read.csv2("02_processed_data/training_test_model.csv")

# pred_mel_eval <- prediction %>% 
#   filter(tribe=="Meliponini") %>%
#   filter(!individual %in% pred_mel_training$individual) %>%
#   group_by(tribe, species) %>%
#   mutate(numb = n()) %>%
#   ungroup()
# 
# write.csv2(pred_mel_eval, "02_processed_data/evaluation_test_model.csv")

pred_mel_eval <- read.csv2("02_processed_data/evaluation_test_model.csv")

### build and predict test model

mel_1 <- lmer(log(proboscisLength) ~ log(intertegularDistance) * genus + (1|numb),
              data = pred_mel_training)

pred_mel_eval$pred <- exp(predict(mel_1, newdata=pred_mel_eval))

pred_mel_eval$diff <- pred_mel_eval$pred - pred_mel_eval$proboscisLength

err_mel <- ggplot(pred_mel_eval) +
  geom_abline(aes(slope=0, intercept=0), linetype="dashed", col="grey40") +
  geom_point(aes(x=intertegularDistance, y=diff), size=2, alpha=0.7) +
  ylim(-2, 2) +
  labs(x="Intertegular Distance [mm]", y="Difference [mm]") +
  theme_bw()

box_mel <- ggplot(pred_mel_eval) +
  geom_violin(aes(x=tribe, y=abs(diff)), fill="grey90") +
  geom_boxplot(aes(x=tribe, y=abs(diff)), fill="lightgreen", width=0.2) +
  ylim(0, 2) +
  labs(x="New model", y="Absolute difference [mm]") +
  theme_bw()

shapiro.test(pred_mel_eval$diff)
t.test(pred_mel_eval$diff, rep(0, 10))

range(pred_mel_eval$diff)
mean(pred_mel_eval$diff)
sd(pred_mel_eval$diff)

quantile(abs(pred_mel_eval$diff), probs=0.90)

### predict Cariveau for evaluation data set

pred_mel_eval$Cariveau <- ITtongue(pred_mel_eval$intertegularDistance, 
                                   pred_mel_eval$family, "tongue")

pred_mel_eval$diff2 <- pred_mel_eval$Cariveau - pred_mel_eval$proboscisLength

err_car <- ggplot(pred_mel_eval) +
  geom_abline(aes(slope=0, intercept=0), linetype="dashed", col="grey40") +
  geom_point(aes(x=intertegularDistance, y=diff2), size=2, alpha=0.7) +
  ylim(-2, 2) +
  labs(x="Intertegular Distance [mm]", y="Difference [mm]") +
  theme_bw()

box_car <- ggplot(pred_mel_eval) +
  geom_violin(aes(x=tribe, y=abs(diff2)), fill="grey90") +
  geom_boxplot(aes(x=tribe, y=abs(diff2)), fill="lightgreen", width=0.2) +
  ylim(0, 2) +
  labs(x="Cariveau et al.", y="Absolute difference [mm]") +
  theme_bw()

shapiro.test(pred_mel_eval$diff2)
t.test(pred_mel_eval$diff2, rep(0, 10))

range(pred_mel_eval$diff2)
mean(pred_mel_eval$diff2)
sd(pred_mel_eval$diff2)

quantile(abs(pred_mel_eval$diff2), probs=0.90)

##### Euglossini excl. Euglossa #####

##### import trait data & create dataframe for model training #####

traits <- read.csv2("02_processed_data/traits.csv")

prediction <- traits %>%
  pivot_wider(id_cols=individual:species, names_from=trait, 
              values_from=value)

##### create linear model for Euglossini #####

prediction <- prediction %>%
  filter(tribe=="Euglossini" & genus != "Euglossa" & genus != "Aglae") 

# pred_eug_training <- prediction %>%
#   sample_n(0.8*132) %>%
#   group_by(tribe, species) %>%
#   mutate(numb = n()) %>%
#   ungroup()
# 
# write.csv2(pred_eug_training, "02_processed_data/training_test_model2.csv")

pred_eug_training <- read.csv2("02_processed_data/training_test_model2.csv")

# pred_eug_eval <- prediction %>%
#   filter(!individual %in% pred_eug_training$individual) %>%
#   group_by(tribe, species) %>%
#   mutate(numb = n()) %>%
#   ungroup()
# 
# write.csv2(pred_eug_eval, "02_processed_data/evaluation_test_model2.csv")

pred_eug_eval <- read.csv2("02_processed_data/evaluation_test_model2.csv")

eug_1 <- lmer(log(proboscisLength) ~ log(intertegularDistance) * genus +  + (1|numb),
              data = pred_eug_training)

fixef(eug_1)

pred <- function(IT, gen){
  
  x <- c()
  
  for(i in 1:length(gen)){
    
    if(gen[i] == "Eufriesea"){
      
      x[i] <- exp(5.31-1.93*log(IT[i]))
      
    } else if(gen[i] == "Eulaema"){
      
      x[i] <- exp(2.23+0.43*log(IT[i]))
      
    } else if(gen[i] == "Exaerete"){
      
      x[i] <- exp(2.12+0.59*log(IT[i]))
      
    } else if(gen[i] == "Euglossa"){
      
      x[i] <- exp(1.79+0.55*log(IT))
      
    }
    
  }
  
  return(x)
  
}

pred_eug_eval$est_newmod <- pred(pred_eug_eval$intertegularDistance, pred_eug_eval$genus)

pred_eug_eval$Cariveau <- ITtongue(pred_eug_eval$intertegularDistance, 
                                   pred_eug_eval$family, "tongue")

pred_eug_eval$diff <- pred_eug_eval$est_newmod - pred_eug_eval$proboscisLength

pred_eug_eval$diff2 <- pred_eug_eval$Cariveau - pred_eug_eval$proboscisLength

shapiro.test(pred_eug_eval$diff)
t.test(pred_eug_eval$diff, rep(0, 10))

range(pred_eug_eval$diff)
mean(pred_eug_eval$diff)
sd(pred_eug_eval$diff)

quantile(abs(pred_eug_eval$diff), probs=0.90)

shapiro.test(pred_eug_eval$diff2)
t.test(pred_eug_eval$diff2, rep(0, 10))

range(pred_eug_eval$diff2)
mean(pred_eug_eval$diff2)
sd(pred_eug_eval$diff2)

quantile(abs(pred_eug_eval$diff2), probs=0.90)

err_eug2 <- ggplot(pred_eug_eval) +
  geom_abline(aes(slope=0, intercept=0), linetype="dashed", col="grey40") +
  geom_point(aes(x=intertegularDistance, y=diff), size=2, alpha=0.7) +
  ylim(-25, 25) +
  labs(x="Intertegular Distance [mm]", y="Difference [mm]") +
  theme_bw()

box_eug2 <- ggplot(pred_eug_eval) +
  geom_violin(aes(x=tribe, y=abs(diff)), fill="grey90") +
  geom_boxplot(aes(x=tribe, y=abs(diff)), fill="lightgreen", width=0.2) +
  ylim(0, 25) +
  labs(x="New model", y="Absolute difference [mm]") +
  theme_bw()

err_car2 <- ggplot(pred_eug_eval) +
  geom_abline(aes(slope=0, intercept=0), linetype="dashed", col="grey40") +
  geom_point(aes(x=intertegularDistance, y=diff2), size=2, alpha=0.7) +
  ylim(-25, 25) +
  labs(x="Intertegular Distance [mm]", y="Difference [mm]") +
  theme_bw()

box_car2 <- ggplot(pred_eug_eval) +
  geom_violin(aes(x=tribe, y=abs(diff2)), fill="grey90") +
  geom_boxplot(aes(x=tribe, y=abs(diff2)), fill="lightgreen", width=0.2) +
  ylim(0, 25) +
  labs(x="Cariveau et al.", y="Absolute difference [mm]") +
  theme_bw()

##### Euglossa #####

##### import trait data & create dataframe for model training #####

traits <- read.csv2("02_processed_data/traits.csv")

prediction <- traits %>%
  pivot_wider(id_cols=individual:species, names_from=trait, 
              values_from=value)

##### create linear model for Euglossini #####

prediction <- prediction %>%
  filter(genus=="Euglossa" & !is.na(subgenus) & !is.na(proboscisLength))

# pred_eug_training <- prediction %>%
#   sample_n(0.8*372) %>%
#   group_by(tribe, species) %>%
#   mutate(numb = n()) %>%
#   ungroup()
# 
# write.csv2(pred_eug_training, "02_processed_data/training_test_model3.csv")

pred_eug_training <- read.csv2("02_processed_data/training_test_model3.csv")

# pred_eug_eval <- prediction %>%
#   filter(!individual %in% pred_eug_training$individual) %>%
#   group_by(tribe, species) %>%
#   mutate(numb = n()) %>%
#   ungroup()
# 
# write.csv2(pred_eug_eval, "02_processed_data/evaluation_test_model3.csv")

pred_eug_eval <- read.csv2("02_processed_data/evaluation_test_model3.csv")

eug_1 <- lmer(log(proboscisLength) ~ log(intertegularDistance) * subgenus +  + (1|numb),
              data = pred_eug_training)

fixef(eug_1)

pred <- function(IT, gen){
  
  x <- c()
  
  for(i in 1:length(gen)){
    
    if(gen[i] == "Euglossa"){
      
      x[i] <- exp(3.32-1.02*log(IT[i]))
      
    } else if(gen[i] == "Glossura"){
      
      x[i] <- exp(0.71+1.85*log(IT[i]))
      
    } else if(gen[i] == "Glossurella"){
      
      x[i] <- exp(2.64+0.07*log(IT[i]))
      
    }
    
  }
  
  return(x)
  
}

pred_eug_eval$est_newmod <- pred(pred_eug_eval$intertegularDistance, pred_eug_eval$subgenus)

pred_eug_eval$Cariveau <- ITtongue(pred_eug_eval$intertegularDistance, 
                                   pred_eug_eval$family, "tongue")

pred_eug_eval$diff <- pred_eug_eval$est_newmod - pred_eug_eval$proboscisLength

pred_eug_eval$diff2 <- pred_eug_eval$Cariveau - pred_eug_eval$proboscisLength

shapiro.test(pred_eug_eval$diff)
t.test(pred_eug_eval$diff, rep(0, 10))

range(pred_eug_eval$diff)
mean(pred_eug_eval$diff)
sd(pred_eug_eval$diff)

quantile(abs(pred_eug_eval$diff), probs=0.90)

shapiro.test(pred_eug_eval$diff2)
t.test(pred_eug_eval$diff2, rep(0, 10))

range(pred_eug_eval$diff2)
mean(pred_eug_eval$diff2)
sd(pred_eug_eval$diff2)

quantile(abs(pred_eug_eval$diff2), probs=0.90)

err_eug3 <- ggplot(pred_eug_eval) +
  geom_abline(aes(slope=0, intercept=0), linetype="dashed", col="grey40") +
  geom_point(aes(x=intertegularDistance, y=diff), size=2, alpha=0.7) +
  ylim(-25, 25) +
  labs(x="Intertegular Distance [mm]", y="Difference [mm]") +
  theme_bw()

box_eug3 <- ggplot(pred_eug_eval) +
  geom_violin(aes(x=genus, y=abs(diff)), fill="grey90") +
  geom_boxplot(aes(x=genus, y=abs(diff)), fill="lightgreen", width=0.2) +
  ylim(0, 25) +
  labs(x="New model", y="Absolute difference [mm]") +
  theme_bw()

err_car3 <- ggplot(pred_eug_eval) +
  geom_abline(aes(slope=0, intercept=0), linetype="dashed", col="grey40") +
  geom_point(aes(x=intertegularDistance, y=diff2), size=2, alpha=0.7) +
  ylim(-25, 25) +
  labs(x="Intertegular Distance [mm]", y="Difference [mm]") +
  theme_bw()

box_car3 <- ggplot(pred_eug_eval) +
  geom_violin(aes(x=genus, y=abs(diff2)), fill="grey90") +
  geom_boxplot(aes(x=genus, y=abs(diff2)), fill="lightgreen", width=0.2) +
  ylim(0, 25) +
  labs(x="Cariveau et al.", y="Absolute difference [mm]") +
  theme_bw()

##### Plot #####
shared_x_label <- textGrob("Intertegular Distance [mm]", gp = gpar(fontsize = 11),
                           hjust=1.5, vjust=0.001)

toprow <- arrangeGrob(
  
  err_mel + labs(title="A", subtitle="New model") + 
    theme(axis.title.x=element_blank()),
  err_car + labs(y="", title="", subtitle="Cariveau et al.") +
    theme(axis.title.x=element_blank()),
  box_mel + labs(x="", title="", subtitle="New model") +
    theme(axis.title.x=element_blank()),
  box_car + labs(x="", y="", title="", subtitle="Cariveau et al.") +
    theme(axis.title.x=element_blank()),
  
  nrow=1
  
)

middlerow <- arrangeGrob(
  
  err_eug2 + labs(title="B", subtitle="New model") +
    theme(axis.title.x=element_blank()),
  err_car2 + labs(y="", title="", subtitle="Cariveau et al.") +
    theme(axis.title.x=element_blank()),
  box_eug2 + labs(x="", title="", subtitle="New model") +
    theme(axis.title.x=element_blank()),
  box_car2 + labs(x="", y="", title="", subtitle="Cariveau et al.") +
    theme(axis.title.x=element_blank()),
  
  nrow=1
  
)

bottomrow <- arrangeGrob(
  
  err_eug3 + labs(title="C", subtitle="New model") +
    theme(axis.title.x=element_blank()),
  err_car3 + labs(y="", title="", subtitle="Cariveau et al.") +
    theme(axis.title.x=element_blank()),
  box_eug3 + labs(x="", title="", subtitle="New model") +
    theme(axis.title.x=element_blank()),
  box_car3 + labs(x="", y="", title="", subtitle="Cariveau et al.") +
    theme(axis.title.x=element_blank()),
  
  nrow=1
  
)

bottomrow_x <- arrangeGrob(
  bottomrow,
  shared_x_label,
  ncol = 1,
  heights = c(10, 1)
)

toprow_x <- arrangeGrob(
  toprow,
  shared_x_label,
  ncol = 1,
  heights = c(10, 1)
)

middlerow_x <- arrangeGrob(
  middlerow,
  shared_x_label,
  ncol = 1,
  heights = c(10, 1)
)

fig2_revised <- arrangeGrob(
  
  toprow_x,
  middlerow_x,
  bottomrow_x
  
)

grid.newpage()
grid.draw(fig2_revised)

ggsave("03_final_data/fig2_revised.png", fig2_revised, width=20, height=20, units="cm")











