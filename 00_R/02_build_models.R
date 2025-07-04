library(tidyverse)
library(BeeIT)
library(MuMIn)
library(gridExtra)
library(lme4)
library(cAIC4)

rm(list=ls())

traits <- read.csv2("02_processed_data/traits.csv")

prediction <- traits %>%
  pivot_wider(id_cols=individual:species, names_from=trait, 
              values_from=value) %>%
  filter(species!="Euglossa sp")

# measurements per species
prediction <- prediction %>%
  group_by(tribe, species) %>%
  mutate(numb = n()) %>%
  ungroup()

range(prediction$numb)

# measured species per tribe
prediction %>%
  group_by(tribe) %>%
  reframe(num = n())

##### predict Cariveau #####

prediction$Cariveau <- ITtongue(prediction$intertegularDistance, 
                                prediction$family, "tongue")

# Meliponini
prediction_mel <- prediction %>% 
  filter(tribe=="Meliponini") %>%
  filter(!is.na(intertegularDistance) & !is.na(proboscisLength))

r_squared_cariveau_mel <- 1 - (sum((prediction_mel$Cariveau 
                                    - prediction_mel$proboscisLength)^2) /
                                 sum((prediction_mel$Cariveau 
                                      - mean(prediction_mel$proboscisLength))^2))

rmse_cariveau_mel <- sqrt(mean((prediction_mel$proboscisLength 
                                - prediction_mel$Cariveau)^2))

range(prediction_mel$Cariveau - prediction_mel$proboscisLength)
mean(prediction_mel$Cariveau - prediction_mel$proboscisLength)
sd(prediction_mel$Cariveau - prediction_mel$proboscisLength)

max(abs(prediction_mel$Cariveau - prediction_mel$proboscisLength)) / 
  prediction_mel$intertegularDistance[max(abs(prediction_mel$Cariveau - 
                                                prediction_mel$proboscisLength))]

# Euglossini
prediction_eug <- prediction %>% 
  filter(tribe=="Euglossini") %>%
  filter(!is.na(intertegularDistance) & !is.na(proboscisLength))

r_squared_cariveau_eug <- 1 - (sum((prediction_eug$Cariveau 
                                    - prediction_eug$proboscisLength)^2) /
                                 sum((prediction_eug$Cariveau 
                                      - mean(prediction_eug$proboscisLength))^2))

rmse_cariveau_eug <- sqrt(mean((prediction_eug$proboscisLength 
                                - prediction_eug$Cariveau)^2))

range(prediction_eug$Cariveau - prediction_eug$proboscisLength)
mean(prediction_eug$Cariveau - prediction_eug$proboscisLength)
sd(prediction_eug$Cariveau - prediction_eug$proboscisLength)

max(abs(prediction_eug$Cariveau - prediction_eug$proboscisLength)) / 
  prediction_eug$intertegularDistance[max(abs(prediction_eug$Cariveau - 
                                                prediction_eug$proboscisLength))]

# Augochlorini
prediction_aug <- prediction %>% 
  filter(tribe=="Augochlorini") %>%
  filter(!is.na(intertegularDistance) & !is.na(proboscisLength))

r_squared_cariveau_aug <- 1 - (sum((prediction_aug$Cariveau 
                                    - prediction_aug$proboscisLength)^2) /
                                 sum((prediction_aug$Cariveau 
                                      - mean(prediction_aug$proboscisLength))^2))

rmse_cariveau_aug <- sqrt(mean((prediction_aug$proboscisLength 
                                - prediction_aug$Cariveau)^2))

range(prediction_aug$Cariveau - prediction_aug$proboscisLength)
mean(prediction_aug$Cariveau - prediction_aug$proboscisLength)
sd(prediction_aug$Cariveau - prediction_aug$proboscisLength)

max(abs(prediction_aug$Cariveau - prediction_aug$proboscisLength)) / 
  prediction_aug$intertegularDistance[max(abs(prediction_aug$Cariveau - 
                                                prediction_aug$proboscisLength))]

##### NEW MODEL FOR MELIPONINI #####

mel_1 <- lmer(log(proboscisLength) ~ log(intertegularDistance) + (1|numb),
            data = prediction_mel)

summary(mel_1)

mel_2 <- lmer(log(proboscisLength) ~ genus + (1|numb), data = prediction_mel)

summary(mel_2)

mel_3 <- lmer(log(proboscisLength) ~ log(intertegularDistance) + genus + (1|numb),
            data = prediction_mel)

summary(mel_3)

mel_4 <- lmer(log(proboscisLength) ~ log(intertegularDistance) * genus + (1|numb),
            data = prediction_mel)

summary(mel_4)

cAIC(mel_1)
cAIC(mel_2)
cAIC(mel_3)
cAIC(mel_4)

r.squaredGLMM(mel_1)
r.squaredGLMM(mel_2)
r.squaredGLMM(mel_3)
r.squaredGLMM(mel_4)

prediction_mel <- prediction_mel %>% 
  mutate(pred_genus=exp(predict(mel_3))) %>%
  mutate(pred_tribe=exp(predict(mel_1)))

fixef(mel_1)
fixef(mel_3)

##### NEW MODEL FOR AUGOCHLORINI #####

aug_1 <- lmer(log(proboscisLength) ~ log(intertegularDistance) + (1|numb),
            data = prediction_aug)

summary(aug_1)

aug_2 <- lmer(log(proboscisLength) ~ genus + (1|numb), data = prediction_aug)

summary(aug_2)

aug_3 <- lmer(log(proboscisLength) ~ log(intertegularDistance) + genus + (1|numb),
            data = prediction_aug)

summary(aug_3)

aug_4 <- lmer(log(proboscisLength) ~ log(intertegularDistance) * genus + (1|numb),
            data = prediction_aug)

summary(aug_4)

cAIC(aug_1)
cAIC(aug_2)
cAIC(aug_3)
cAIC(aug_4)

r.squaredGLMM(aug_1)
r.squaredGLMM(aug_2)
r.squaredGLMM(aug_3)
r.squaredGLMM(aug_4)

prediction_aug <- prediction_aug %>% 
  mutate(pred_genus=exp(predict(aug_3))) %>%
  mutate(pred_tribe=exp(predict(aug_1)))

fixef(aug_1)
fixef(aug_3)

##### NEW MODELS FOR ALL EUGLOSSINI #####

eug_1 <- lmer(log(proboscisLength) ~ log(intertegularDistance) + (1|numb),
            data = prediction_eug)

summary(eug_1)

eug_2 <- lmer(log(proboscisLength) ~ genus + (1|numb), data = prediction_eug)

summary(eug_2)

eug_3 <- lmer(log(proboscisLength) ~ log(intertegularDistance) + genus + (1|numb),
            data = prediction_eug)

summary(eug_3)

eug_4 <- lmer(log(proboscisLength) ~ log(intertegularDistance) * genus + (1|numb),
            data = prediction_eug)

summary(eug_4)

cAIC(eug_1)
cAIC(eug_2)
cAIC(eug_3)
cAIC(eug_4)

r.squaredGLMM(eug_1)
r.squaredGLMM(eug_2)
r.squaredGLMM(eug_3)
r.squaredGLMM(eug_4)

prediction_eug$pred_genus <- exp(predict(eug_4))
prediction_eug$pred_tribe <- exp(predict(eug_1))

fixef(eug_1)
fixef(eug_3)

##### NEW MODELS FOR EUGLOSSINI WITHOUT EUGLOSSA #####

prediction_eug.without_euglossa <- prediction_eug %>%
  filter(genus!="Euglossa" & genus!="Aglae") %>%
  filter(!is.na(intertegularDistance) & !is.na(proboscisLength))

eug_1.without_euglossa <- lmer(log(proboscisLength) ~ log(intertegularDistance) + (1|numb),
            data = prediction_eug.without_euglossa)

summary(eug_1.without_euglossa)

eug_2.without_euglossa <- lmer(log(proboscisLength) ~ genus + (1|numb), data = prediction_eug.without_euglossa)

summary(eug_2.without_euglossa)

eug_3.without_euglossa <- lmer(log(proboscisLength) ~ log(intertegularDistance) + genus + (1|numb),
            data = prediction_eug.without_euglossa)

summary(eug_3.without_euglossa)

eug_4.without_euglossa <- lmer(log(proboscisLength) ~ log(intertegularDistance) * genus + (1|numb),
            data = prediction_eug.without_euglossa)

summary(eug_4.without_euglossa)

cAIC(eug_1.without_euglossa)
cAIC(eug_2.without_euglossa)
cAIC(eug_3.without_euglossa)
cAIC(eug_4.without_euglossa)

r.squaredGLMM(eug_1.without_euglossa)
r.squaredGLMM(eug_2.without_euglossa)
r.squaredGLMM(eug_3.without_euglossa)
r.squaredGLMM(eug_4.without_euglossa)

prediction_eug.without_euglossa$pred_genus <- exp(predict(eug_4.without_euglossa))
prediction_eug.without_euglossa$pred_tribe <- exp(predict(eug_1.without_euglossa))

fixef(eug_1.without_euglossa)
fixef(eug_4.without_euglossa)

##### NEW MODELS FOR EUGLOSSA INCLUDING SUBGENERA #####

prediction_eug.only_euglossa <- prediction_eug %>%
  filter(genus=="Euglossa") %>%
  filter(!is.na(intertegularDistance) & !is.na(proboscisLength) & !is.na(subgenus))

eug_1.only_euglossa <- lmer(log(proboscisLength) ~ log(intertegularDistance) + (1|numb),
            data = prediction_eug.only_euglossa)

summary(eug_1.only_euglossa)

eug_2.only_euglossa <- lmer(log(proboscisLength) ~ subgenus + (1|numb), data = prediction_eug.only_euglossa)

summary(eug_2.only_euglossa)

eug_3.only_euglossa <- lmer(log(proboscisLength) ~ log(intertegularDistance) + subgenus + (1|numb),
            data = prediction_eug.only_euglossa)

summary(eug_3.only_euglossa)

eug_4.only_euglossa <- lmer(log(proboscisLength) ~ log(intertegularDistance) * subgenus + (1|numb),
            data = prediction_eug.only_euglossa)

summary(eug_4.only_euglossa)

cAIC(eug_1.only_euglossa)
cAIC(eug_2.only_euglossa)
cAIC(eug_3.only_euglossa)
cAIC(eug_4.only_euglossa)

r.squaredGLMM(eug_1.only_euglossa)
r.squaredGLMM(eug_2.only_euglossa)
r.squaredGLMM(eug_3.only_euglossa)
r.squaredGLMM(eug_4.only_euglossa)

prediction_eug.only_euglossa$pred_genus <- exp(predict(eug_4.only_euglossa))
prediction_eug.only_euglossa$pred_tribe <- exp(predict(eug_1.only_euglossa))

fixef(eug_1.only_euglossa)
fixef(eug_4.only_euglossa)

##### PLOTS #####

prediction <- rbind(prediction_mel, prediction_aug, prediction_eug.without_euglossa, prediction_eug.only_euglossa) %>%
  mutate(tribe_plot=ifelse(genus=="Euglossa", "Euglossini (only Euglossa)", tribe)) %>%
  mutate(genus_plot=ifelse(tribe_plot=="Euglossini (only Euglossa)", subgenus, genus))

facet_titles <- as_labeller(
  c("Augochlorini" = paste0("Augochlorini; R² = ",
                            round(r.squaredGLMM(aug_1)[2], 2)),
    "Meliponini" = paste0("Meliponini; R² = ",
                          round(r.squaredGLMM(mel_1)[2], 2)),
    "Euglossini" = paste0("Euglossini; R² = ",
                          round(r.squaredGLMM(eug_1)[2], 2)))
)

plot_tribe <- ggplot(prediction, aes(x=intertegularDistance)) +
  scale_color_manual(breaks=c("Augochlorini", "Meliponini", "Euglossini"),
                     values=c("magenta4", "turquoise4", "darkgreen")) +
  scale_fill_manual(breaks=c("Augochlorini", "Meliponini", "Euglossini"),
                    values=c("magenta4", "turquoise4", "darkgreen")) +
  geom_point(aes(y=proboscisLength, col=tribe), alpha=0.06, size=2) +
  geom_smooth(data=prediction, 
              aes(y=pred_tribe, col=tribe, fill=tribe), method="lm",
              show.legend=F, se=F) +
  geom_smooth(aes(y=Cariveau), col="black", linetype="dashed") +
  theme_bw() +
  facet_wrap(.~factor(tribe, 
                      levels=c("Augochlorini", "Meliponini", "Euglossini")),
             scales="free",
             labeller=facet_titles) +
  labs(x="Intertegular Distance [mm]", y="Proboscis Length [mm]", 
       col="Tribe", fill="Tribe") +
  theme(strip.background=element_rect(fill="grey95"),
        legend.position="none")

facet_titles2 <- as_labeller(
  c("Augochlorini" = paste0("Augochlorini; R² = ",
                            round(r.squaredGLMM(aug_3)[2], 2)),
    "Meliponini" = paste0("Meliponini; R² = ",
                          round(r.squaredGLMM(mel_3)[2], 2)),
    "Euglossini" = paste0("Euglossini; R² = ",
                          round(r.squaredGLMM(eug_4.without_euglossa)[2], 2),
                          " / (Euglossa; R² = ",
                          round(r.squaredGLMM(eug_4.only_euglossa)[2], 2), ")"))
)

plot_genus <- ggplot(prediction, aes(x=intertegularDistance)) +
  scale_color_manual(breaks=c("Augochlorini", "Meliponini", "Euglossini"),
                     values=c("magenta4", "turquoise4", "darkgreen")) +
  geom_point(aes(y=proboscisLength, col=tribe), alpha=0.06, size=2) +
  geom_smooth(aes(y=pred_genus, fill=genus_plot, col=tribe), method="lm", se=F) +
  theme_bw() +
  facet_wrap(.~factor(tribe, 
                      levels=c("Augochlorini", "Meliponini", "Euglossini")),
             scales="free",
             labeller=facet_titles2) +
  labs(x="Intertegular Distance [mm]", y="Proboscis Length [mm]") +
  theme(strip.background=element_rect(fill="grey95"),
        legend.position="none")

fig1 <- grid.arrange(
  plot_tribe + labs(title="A"), 
  plot_genus + labs(title="B")
)

ggsave("03_final_data/fig_1.png", fig1, dpi=600)
