library(naniar)

miss_var_summary(PFAS_ALDUST_DEMO_2005) %>% data.frame()
miss_var_table(PFAS_ALDUST_DEMO_2005)

### deal with missing values
full <- PFAS_ALDUST_DEMO_2005 %>% left_join(BIOPRO_D.1, by = "SEQN") %>% 
  select(SDMVPSU, SDMVSTRA, WTSA2YR, 
         LBXPFOA, LBXPFOS, LBXPFHS, LBXMPAH, LBXPFDE, LBXPFNA, 
         RIDAGEYR, gender, race, edu, INDFMPIR, 
         eatout.cat, shellfish, fish, tapwater, BMXBMI,
         carpet.cat4, country, military, 
         period, pregnant, pregnant.times, live.birth, breastfed, children.breastfed,
         LBXSCR) %>% 
  mutate(Standard.Creatinine = -0.016 + 0.978*LBXSCR)


full.miss <- full %>% select(LBXPFOA, LBXPFOS, LBXPFHS, LBXMPAH, LBXPFDE, LBXPFNA, 
                             race, gender, RIDAGEYR, edu, INDFMPIR, 
                             eatout.cat, shellfish, fish, tapwater, BMXBMI, 
                             carpet.cat4, country, military,
                             period, pregnant.times, children.breastfed,
                             Standard.Creatinine) %>%
  dplyr::rename(PFOA = LBXPFOA, PFOS = LBXPFOS, PFHxS = LBXPFHS, MeFOSAA = LBXMPAH, PFDA = LBXPFDE, PFNA = LBXPFNA,
                `Race/ethnicity` = race, Gender = gender, Age = RIDAGEYR, Education = edu, `Family PIR` = INDFMPIR, 
                `Eating out per week` = eatout.cat, 
                `Eating shellfish during past 30 days` = shellfish, 
                `Eating fish during past 30 days` = fish, 
                `Tap water source` = tapwater, BMI = BMXBMI, 
                `Type of floor covering` = carpet.cat4, `Country of birth` = country, `Veteran/military status` = military,
                `Had at least one period in the past 12 months` = period, 
                `Number of pregnancies` = pregnant.times,  
                `Number of children breastfed for at least 1 month` = children.breastfed,
                `Standard Creatinine` = Standard.Creatinine)

workdir <- " "
par(mai = c(1.5, 1.5, 1.5, 1.5))
tiff(file.path(workdir, "MissingData.tiff"), width = 12, height = 8, units = 'in', res = 600)
full.miss %>% vis_miss(cluster = TRUE, sort_miss = TRUE)
dev.off()
full.miss %>% is.na() %>% colSums()

full.miss %>% arrange(PFOA) %>% vis_miss()

full %>% vis_miss(sort_miss = TRUE)

full %>% gg_miss_case(facet = carpet.cat4)
full %>% gg_miss_case(facet = edu)

tiff(file.path(workdir, "MissingData1.tiff"), width = 12, height = 8, units = 'in', res = 600)
full.miss %>% gg_miss_var()
dev.off()


full.miss %>% rename(carpet.cat4 = `Type of floor covering`) %>% gg_miss_var(facet = carpet.cat4)


gg_miss_upset(full)
gg_miss_fct(x = full, fct = carpet.cat4)


tiff(file.path(workdir, "MissingData2.tiff"), width = 12, height = 8, units = 'in', res = 600)
gg_miss_fct(x = full.miss %>% rename(Race = `Race/ethnicity`), fct = Race)
dev.off()


tiff(file.path(workdir, "MissingData3.tiff"), width = 12, height = 8, units = 'in', res = 600)
gg_miss_fct(x = full.miss, fct = Education)
dev.off()

gg_miss_fct(x = full.miss %>% rename(country = `Country of birth`), fct = country)

as_shadow(full) %>% View()
bind_shadow(full, only_miss = TRUE) %>% View()

bind_shadow(full, only_miss = TRUE) %>% group_by(carpet.cat4) %>%
  summarise(PFOA_mean = mean(LBXPFOA, na.rm = TRUE),
            PFOS_mean = mean(LBXPFOS, na.rm = TRUE),
            PFHS_mean = mean(LBXPFHS, na.rm = TRUE),
            MPAH_mean = mean(LBXMPAH, na.rm = TRUE),
            PFDE_mean = mean(LBXPFDE, na.rm = TRUE),
            PFNA_mean = mean(LBXPFNA, na.rm = TRUE))

library(ggpubr)

fig1 <- bind_shadow(full, only_miss = TRUE) %>% ggplot(aes(x = LBXPFOA, color = carpet.cat4_NA)) + geom_density() + 
  labs(color = "Type of floor covering") + xlab("PFOA Concentrations (ng/ml)")
fig2 <- bind_shadow(full, only_miss = TRUE) %>% ggplot(aes(x = LBXPFOS, color = carpet.cat4_NA)) + geom_density() +
  labs(color = "Type of floor covering") + xlab("PFOS Concentrations (ng/ml)")
fig3 <- bind_shadow(full, only_miss = TRUE) %>% ggplot(aes(x = LBXPFHS, color = carpet.cat4_NA)) + geom_density() +
  labs(color = "Type of floor covering") + xlab("PFHxS Concentrations (ng/ml)")
fig4 <- bind_shadow(full, only_miss = TRUE) %>% ggplot(aes(x = LBXMPAH, color = carpet.cat4_NA)) + geom_density() +
  labs(color = "Type of floor covering") + xlab("MeFOSAA Concentrations (ng/ml)")
fig5 <- bind_shadow(full, only_miss = TRUE) %>% ggplot(aes(x = LBXPFDE, color = carpet.cat4_NA)) + geom_density() +
  labs(color = "Type of floor covering") + xlab("PFDA Concentrations (ng/ml)")
fig6 <- bind_shadow(full, only_miss = TRUE) %>% ggplot(aes(x = LBXPFNA, color = carpet.cat4_NA)) + geom_density() +
  labs(color = "Type of floor covering") + xlab("PFNA Concentrations (ng/ml)")


tiff(file.path(workdir, "MissingData4.tiff"), width = 24, height = 8, units = 'in', res = 600)
fig <- ggarrange(fig1, fig2, fig3, fig4, fig5, fig6,
                 labels = c("A", "B", "C", "D", "E", "F"),
                 ncol = 3, nrow = 2)
fig
dev.off()

bind_shadow(full, only_miss = TRUE) %>% ggplot(aes(x = LBXPFOA, color = carpet.cat4_NA)) + geom_density() 

bind_shadow(full, only_miss = TRUE) %>% ggplot(aes(x = LBXPFOA, color = carpet.cat4_NA)) + geom_density() +
  facet_wrap(~race)
bind_shadow(full, only_miss = TRUE) %>% ggplot(aes(x = carpet.cat4_NA, y = LBXPFOA)) + geom_boxplot()
bind_shadow(full, only_miss = TRUE) %>% ggplot(aes(x = carpet.cat4_NA, y = LBXPFOA)) + geom_boxplot() +
  facet_wrap(~RIDRETH1)

ggplot(full, aes(x = carpet.cat4, y = LBXPFOA)) + geom_miss_point()
ggplot(full, aes(x = carpet.cat4, y = LBXPFOA)) + geom_miss_point() + facet_wrap(~edu)
bind_shadow(full) %>% ggplot(aes(x = carpet.cat4, y = LBXPFOA)) + geom_miss_point() + facet_wrap(~military_NA)
bind_shadow(full) %>% ggplot(aes(x = carpet.cat4, y = LBXPFOA)) + geom_miss_point() + facet_wrap(edu~military_NA)

dev.set(dev.next())   ### start showing plot in R


full.m <- full %>% mutate(missing_carpet = is.na(carpet.cat4))

missing_carpet_female <- full.m %>% filter(gender == "female") %>% pull(missing_carpet)
missing_carpet_male <- full.m %>% filter(gender == "male") %>% pull(missing_carpet)
t.test(missing_carpet_female, missing_carpet_male)  # not missing at random with respect to gender

missing_carpet_white <- full.m %>% filter(race == "White") %>% pull(missing_carpet)
missing_carpet_black <- full.m %>% filter(race == "Black") %>% pull(missing_carpet)
missing_carpet_hispanic <- full.m %>% filter(race == "Hispanic") %>% pull(missing_carpet)
missing_carpet_others <- full.m %>% filter(race == "Other Race") %>% pull(missing_carpet)

library(VIM)
par(mai = c(2.5, 1, 0.42, 0.42))
full %>% select(race, carpet.cat4) %>% rename(`type of floor covering` = carpet.cat4) %>% spineMiss(xlab = "")
full %>% select(edu, carpet.cat4) %>% rename(`type of floor covering` = carpet.cat4) %>% spineMiss(xlab = "")
full %>% select(country, carpet.cat4) %>% rename(`type of floor covering` = carpet.cat4) %>% spineMiss(xlab = "")
full %>% select(military, carpet.cat4) %>% rename(`type of floor covering` = carpet.cat4) %>% spineMiss(xlab = "")
full %>% select(country, military) %>% rename(`veteran/military status` = military) %>% spineMiss(xlab = "")
full %>% select(gender, military) %>% rename(`veteran/military status` = military) %>% spineMiss(xlab = "")
full %>% select(race, military) %>% rename(`veteran/military status` = military) %>% spineMiss(xlab = "")
full %>% select(edu, military) %>% rename(`veteran/military status` = military) %>% spineMiss(xlab = "")
