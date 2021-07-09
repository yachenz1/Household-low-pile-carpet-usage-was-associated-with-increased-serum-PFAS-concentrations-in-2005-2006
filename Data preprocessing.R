library(SASxport)
library(Hmisc)
library(dplyr)
library(naniar)
library(mice)
library(mitools)
library(survey)
library(broom)

############################################# NHANES 2005-2006 ###############################################

# Allergens - Household Dust
ALDUST_D <- read.xport("ALDUST_D.XPT") 
ALDUST_D.1 <- ALDUST_D[, c("SEQN", "AAXFLTYP", "LBXDWT")]


# Polyfluoroalkyl Chemicals
PFC_D <- read.xport("PFC_D.XPT")


# Body Measures
BMX_D <- read.xport("BMX_D.XPT")
BMXKeepVars_D <- c("SEQN", "BMXBMI")
BMX_D.1 <- BMX_D[, BMXKeepVars_D]


# Reproductive Health Questionnaire
RHQ_D <- read.xport("RHQ_D.XPT")
RHQKeepVars_D <- c("SEQN",
                   "RHQ031", # at least one period in past 12 months
                   "RHQ131", # Ever been pregnant
                   "RHQ160", # how many times have been pregnant
                   "RHQ166", # how many vaginal deliveries
                   "RHQ169", # how many cesarean deliveries
                   "RHQ171", # how many deliveries live birth result
                   "RHQ210", # breastfed any of your children
                   "RHD230" # of children breastfed at least one month
)
RHQ_D.1 <- RHQ_D[, RHQKeepVars_D]


# Dietary Interview - Total Nutrient Intakes, First Day
DR1TOT_D <- read.xport("DR1TOT_D.XPT")
DR1TOTKeepVars_D <- c("SEQN",
                      "DR1DRSTZ",
                      "DR1TWS", # tap water source
                      "DR1.330Z", # tap water drank yesterday
                      "DRD340", # shellfish eaten in the past 30 days
                      "DRD360") # fish eaten in the past 30 days
DR1TOT_D.1 <- DR1TOT_D[, DR1TOTKeepVars_D]


# Dietary Behavior Questionnaire
DBQ_D <- read.xport("DBQ_D.XPT")
DBQKeepVars_D <- c("SEQN", 
                   "DBD091") # number of meals eaten outside
DBQ_D.1 <- DBQ_D[, DBQKeepVars_D]


# Standard Biochemistry Profile
BIOPRO_D <- read.xport("BIOPRO_D.XPT")
BIOPRO_D.1 <- BIOPRO_D[, c("SEQN", "LBXSCR")] 


PFAS_ALDUST_2005 <- merge(PFC_D, ALDUST_D.1, by = "SEQN", all = F)


# Demographic Variables & Sample Weights
DEMO_D <- read.xport("DEMO_D.XPT")
DemoKeepVars_D <- c("SEQN",  # unique person identifier (merge variable)
                    "SDMVPSU", # primary sampling unit variable, used in complex design
                    "SDMVSTRA", # strata variable
                    "RIDRETH1", # person race
                    "RIDAGEYR", # Age in years at screening
                    "RIAGENDR", # gender
                    "DMQMILIT", # Served active duty in US Armed Forces
                    "DMDBORN", # Country of birth
                    "DMDEDUC3", # Education level - Children/Youth 6-19
                    "DMDEDUC2", # Education level - Adults 20+
                    "INDFMPIR" # Ratio of family income to poverty
)
DEMO_D.1 <- DEMO_D[, DemoKeepVars_D]


PFAS_ALDUST_DEMO_2005 <- merge(PFAS_ALDUST_2005, DEMO_D.1, by = "SEQN", all = F) %>%
  left_join(DR1TOT_D.1, by = "SEQN") %>% left_join(BMX_D.1, by = "SEQN") %>%
  left_join(RHQ_D.1, by = "SEQN") %>% left_join(DBQ_D.1, by = "SEQN") 


# Assign race/ethnicity levels
race.lev <- c("Hispanic", "White", "Black", "Other Race")
PFAS_ALDUST_DEMO_2005$race <- case_when(PFAS_ALDUST_DEMO_2005$RIDRETH1 == 1 | PFAS_ALDUST_DEMO_2005$RIDRETH1 == 2 ~ race.lev[1], 
                                        PFAS_ALDUST_DEMO_2005$RIDRETH1 == 3 ~ race.lev[2],
                                        PFAS_ALDUST_DEMO_2005$RIDRETH1 == 4 ~ race.lev[3],
                                        PFAS_ALDUST_DEMO_2005$RIDRETH1 == 5 ~ race.lev[4])
PFAS_ALDUST_DEMO_2005$race <- relevel(factor(PFAS_ALDUST_DEMO_2005$race), ref = 'Hispanic')
table(PFAS_ALDUST_DEMO_2005$RIDRETH1)


# Assign education levels
DMDEDUC3 <- PFAS_ALDUST_DEMO_2005$DMDEDUC3
DMDEDUC2 <- PFAS_ALDUST_DEMO_2005$DMDEDUC2

edulev <- c("Less than college", "Some college", "College graduate or above")
PFAS_ALDUST_DEMO_2005$edu <- case_when(DMDEDUC3 <= 14 | DMDEDUC3 == 66 | DMDEDUC2 %in% 1:3 ~ edulev[1], 
                                       DMDEDUC3 == 15 | DMDEDUC2 == 4 ~ edulev[2], 
                                       DMDEDUC2 == 5 ~ edulev[3])
PFAS_ALDUST_DEMO_2005$edu <- relevel(factor(PFAS_ALDUST_DEMO_2005$edu), ref = 'Less than college')
table(PFAS_ALDUST_DEMO_2005$edu)
PFAS_ALDUST_DEMO_2005 %>% filter(is.na(edu)) %>% nrow()


# Assign gender levels
gendlev <- c("male", "female")
PFAS_ALDUST_DEMO_2005$gender <- case_when(PFAS_ALDUST_DEMO_2005$RIAGENDR == 1 ~ gendlev[1], 
                                          PFAS_ALDUST_DEMO_2005$RIAGENDR == 2 ~ gendlev[2])
PFAS_ALDUST_DEMO_2005$gender <- relevel(factor(PFAS_ALDUST_DEMO_2005$gender), ref = 'female')
table(PFAS_ALDUST_DEMO_2005$RIAGENDR)
#PFAS_ALDUST_DEMO_2005.1$RIAGENDR <- factor(PFAS_ALDUST_DEMO_2005.1$RIAGENDR, levels = 1:2, labels = c("male", "female"))


# Assign period levels
periodlev <- c("Yes", "No")
PFAS_ALDUST_DEMO_2005$period <- case_when(PFAS_ALDUST_DEMO_2005$RHQ031 == 1 ~ periodlev[1],
                                          PFAS_ALDUST_DEMO_2005$RHQ031 == 2 | PFAS_ALDUST_DEMO_2005$gender == "male" | PFAS_ALDUST_DEMO_2005$RIDAGEYR >= 55 ~ periodlev[2])
PFAS_ALDUST_DEMO_2005$period <- relevel(factor(PFAS_ALDUST_DEMO_2005$period), ref = "No")


# Assign pregnant levels
pregnantlev <- c("Yes", "No")
PFAS_ALDUST_DEMO_2005$pregnant <- case_when(PFAS_ALDUST_DEMO_2005$RHQ131 == 1 ~ pregnantlev[1],
                                            PFAS_ALDUST_DEMO_2005$RHQ131 == 2 | PFAS_ALDUST_DEMO_2005$gender == "male" ~ pregnantlev[2])
PFAS_ALDUST_DEMO_2005$pregnant <- relevel(factor(PFAS_ALDUST_DEMO_2005$pregnant), ref = "No")

# Pregnant times
PFAS_ALDUST_DEMO_2005$pregnant.times <- ifelse(PFAS_ALDUST_DEMO_2005$pregnant == "No", 0, PFAS_ALDUST_DEMO_2005$RHQ160)

PFAS_ALDUST_DEMO_2005 <- replace_with_na(PFAS_ALDUST_DEMO_2005, replace = list(pregnant.times = 99))
PFAS_ALDUST_DEMO_2005 %>% select(RHQ160, pregnant, pregnant.times) %>% View()

# Live births
PFAS_ALDUST_DEMO_2005$live.birth <- ifelse(PFAS_ALDUST_DEMO_2005$pregnant == "No", 0, PFAS_ALDUST_DEMO_2005$RHQ171)
PFAS_ALDUST_DEMO_2005 %>% select(pregnant, live.birth, RHQ171) %>% View()

# Breastfeeding
PFAS_ALDUST_DEMO_2005$breastfed <- case_when(PFAS_ALDUST_DEMO_2005$RHQ210 == 1 ~ "Yes",
                                             PFAS_ALDUST_DEMO_2005$RHQ210 == 2 | PFAS_ALDUST_DEMO_2005$pregnant == "No" ~ "No")
PFAS_ALDUST_DEMO_2005$breastfed <- relevel(factor(PFAS_ALDUST_DEMO_2005$breastfed), ref = "No")
PFAS_ALDUST_DEMO_2005 %>% select(RIAGENDR, RHQ131, pregnant, breastfed, RHQ210, RIAGENDR) %>% View()


PFAS_ALDUST_DEMO_2005 <- PFAS_ALDUST_DEMO_2005 %>% mutate(children.breastfed = RHD230) %>% 
  mutate(children.breastfed = replace(children.breastfed, breastfed == "No", 0)) %>%
  mutate(breastfed = replace(breastfed, is.na(breastfed) & children.breastfed > 0, "Yes"))


PFAS_ALDUST_DEMO_2005 %>% select(RIAGENDR, RHQ131, pregnant, breastfed, RHQ210, children.breastfed, RHD230) %>% View()



TapWaterSource <- c("Community Supply", "Don't drink tap water", "Other")
PFAS_ALDUST_DEMO_2005$tapwater <- case_when(PFAS_ALDUST_DEMO_2005$DR1TWS == 1 ~ TapWaterSource[1],
                                            PFAS_ALDUST_DEMO_2005$DR1TWS == 4 ~ TapWaterSource[2],
                                            PFAS_ALDUST_DEMO_2005$DR1TWS %in% c(2, 3, 91) ~ TapWaterSource[3])
PFAS_ALDUST_DEMO_2005$tapwater <- relevel(factor(PFAS_ALDUST_DEMO_2005$tapwater), ref = "Don't drink tap water")


militarylev <- c("Yes", "No")
PFAS_ALDUST_DEMO_2005$military <- case_when( PFAS_ALDUST_DEMO_2005$DMQMILIT == 1 ~ militarylev[1],
                                             PFAS_ALDUST_DEMO_2005$DMQMILIT == 2 ~ militarylev[2])
PFAS_ALDUST_DEMO_2005$military <- relevel(factor(PFAS_ALDUST_DEMO_2005$military), ref = "No")

countrylev <- c("US", "Foreign")
PFAS_ALDUST_DEMO_2005$country <- case_when(PFAS_ALDUST_DEMO_2005$DMDBORN == 1 ~ countrylev[1],
                                           PFAS_ALDUST_DEMO_2005$DMDBORN %in% 2:3 ~ countrylev[2])
PFAS_ALDUST_DEMO_2005$country <- relevel(factor(PFAS_ALDUST_DEMO_2005$country), ref = "Foreign")

eatoutlev <- c("No", "Yes")
PFAS_ALDUST_DEMO_2005$eatout.cat <- case_when(PFAS_ALDUST_DEMO_2005$DBD091 == 0 ~ "No",
                                              PFAS_ALDUST_DEMO_2005$DBD091 %in% c(1:21, 5555, 6666) ~ "Yes")
PFAS_ALDUST_DEMO_2005$eatout.cat <- relevel(factor(PFAS_ALDUST_DEMO_2005$eatout.cat), ref = 'No')
table(PFAS_ALDUST_DEMO_2005$eatout.cat)
PFAS_ALDUST_DEMO_2005 %>% filter(is.na(eatout.cat)) %>% nrow()


shellfish.lev <- c("Yes", "No")
PFAS_ALDUST_DEMO_2005$shellfish <- case_when( PFAS_ALDUST_DEMO_2005$DRD340 == 1 ~ "Yes",
                                              PFAS_ALDUST_DEMO_2005$DRD340 == 2 ~ "No")
PFAS_ALDUST_DEMO_2005$shellfish <- relevel(factor(PFAS_ALDUST_DEMO_2005$shellfish), ref = 'No')
table(PFAS_ALDUST_DEMO_2005$DRD340)
label(PFAS_ALDUST_DEMO_2005$DRD340) <- "shellfish eaten in the past 30 days"

fish.lev <- c("Yes", "No")
PFAS_ALDUST_DEMO_2005$fish <- case_when(PFAS_ALDUST_DEMO_2005$DRD360 == 1 ~ "Yes",
                                        PFAS_ALDUST_DEMO_2005$DRD360 == 2 ~ "No")
PFAS_ALDUST_DEMO_2005$fish <- relevel(factor(PFAS_ALDUST_DEMO_2005$fish), ref = 'No')
table(PFAS_ALDUST_DEMO_2005$DRD360)
label(PFAS_ALDUST_DEMO_2005$DRD360) <- "fish eaten in the past 30 days"


# Keep the four categories in Carpeting
carpet <- PFAS_ALDUST_DEMO_2005$AAXFLTYP
carpet.lev4 <- c("Low pile carpet or rug", "Medium to high pile carpet or rug", "Smooth surface", 
                 "Combination of carpet and smooth surface")
PFAS_ALDUST_DEMO_2005$carpet.cat4 <- case_when(carpet == 1 ~ carpet.lev4[1],
                                               carpet == 2 ~ carpet.lev4[2],
                                               carpet == 3 ~ carpet.lev4[3],
                                               carpet == 4 ~ carpet.lev4[4])
PFAS_ALDUST_DEMO_2005$carpet.cat4 <- relevel(factor(PFAS_ALDUST_DEMO_2005$carpet.cat4), ref = 'Smooth surface')
table(PFAS_ALDUST_DEMO_2005$carpet.cat4)

PFAS_ALDUST_DEMO_2005 %>% filter(!is.na(carpet.cat4)) %>% nrow() # 1698

PFAS_ALDUST_DEMO_2005 %>% filter(!is.na(RIAGENDR) & !is.na(RIDAGEYR) & !is.na(RIDRETH1)) %>% 
  filter(!is.na(edu) & !is.na(INDFMPIR) & !is.na(carpet.cat4) & !is.na(eatout.cat)) %>%
  filter(!is.na(DRD340) & !is.na(DRD360)) %>% nrow() 

PFAS_ALDUST_DEMO_2005 %>% filter(!is.na(RIAGENDR) & !is.na(RIDAGEYR) & !is.na(RIDRETH1)) %>% 
  filter(!is.na(edu) & !is.na(INDFMPIR) & !is.na(carpet.cat4) & !is.na(eatout.cat)) %>%
  filter(!is.na(DRD340) & !is.na(DRD360) & !is.na(BMXBMI)) %>% nrow() 

