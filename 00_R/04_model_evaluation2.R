library(tidyverse)
library(BeeIT)
library(MuMIn)
library(gridExtra)
library(lme4)

rm(list=ls())

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

err_eug <- ggplot(pred_eug_eval) +
  geom_abline(aes(slope=0, intercept=0), linetype="dashed", col="grey40") +
  geom_point(aes(x=intertegularDistance, y=diff), size=2, alpha=0.7) +
  ylim(-25, 25) +
  labs(x="Intertegular Distance [mm]", y="Difference [mm]") +
  theme_bw()

box_eug <- ggplot(pred_eug_eval) +
  geom_violin(aes(x=tribe, y=abs(diff)), fill="grey90") +
  geom_boxplot(aes(x=tribe, y=abs(diff)), fill="lightgreen", width=0.2) +
  ylim(0, 25) +
  labs(x="New model", y="Absolute difference [mm]") +
  theme_bw()

err_car <- ggplot(pred_eug_eval) +
  geom_abline(aes(slope=0, intercept=0), linetype="dashed", col="grey40") +
  geom_point(aes(x=intertegularDistance, y=diff2), size=2, alpha=0.7) +
  ylim(-25, 25) +
  labs(x="Intertegular Distance [mm]", y="Difference [mm]") +
  theme_bw()

box_car <- ggplot(pred_eug_eval) +
  geom_violin(aes(x=tribe, y=abs(diff2)), fill="grey90") +
  geom_boxplot(aes(x=tribe, y=abs(diff2)), fill="lightgreen", width=0.2) +
  ylim(0, 25) +
  labs(x="Cariveau et al.", y="Absolute difference [mm]") +
  theme_bw()

fig3 <- grid.arrange(
  err_eug + labs(title="New model"),
  err_car + labs(title="Cariveau et al."),
  box_eug + labs(x=""),
  box_car + labs(x=""),
  ncol=2
)

ggsave("03_final_data/fig_S2.png", fig3, dpi=600)



