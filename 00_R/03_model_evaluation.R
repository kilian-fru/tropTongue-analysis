library(tidyverse)
library(BeeIT)
library(MuMIn)
library(gridExtra)

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

fig2 <- grid.arrange(
  err_mel + labs(title="New model"),
  err_car + labs(title="Cariveau et al."),
  box_mel + labs(x=""),
  box_car + labs(x=""),
  ncol=2
)

ggsave("03_final_data/fig_2.png", fig2, dpi=600)
