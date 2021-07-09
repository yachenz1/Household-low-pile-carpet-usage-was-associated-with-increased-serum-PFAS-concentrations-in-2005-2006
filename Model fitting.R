### Multiple imputation with k=100

k = 100; p <- 22 + 1  # change p to 22 when run the analysis for adolescents only

# remove the irrelevant variables before the imputation
vars_by_NAs <- full %>% select(-c(SDMVPSU, SDMVSTRA, WTSA2YR, pregnant, live.birth, breastfed, LBXSCR)) %>% 
  is.na() %>% colSums() %>% sort(decreasing = FALSE) %>% names() # sort by missing rate (increasing)
multiimp <- mice(full %>% select(vars_by_NAs), m = k, seed = 1)
betas1 <- ses1 <- betas2 <- ses2 <- betas3 <- ses3 <- betas4 <- ses4 <- betas5 <- ses5 <- betas6 <- ses6 <- matrix(NA, k, p) # stats for adjusted models
betas1.c <- ses1.c <- betas2.c <- ses2.c <- betas3.c <- ses3.c <- betas4.c <- ses4.c <- betas5.c <- ses5.c <- betas6.c <- ses6.c <- matrix(NA, k, 4) # stats for crude models
cov.PFOA <- cov.PFOS <- cov.PFHxS <- cov.MeFOSAA <- cov.PFDA <- cov.PFNA <- array(dim = c(p,p,k))


stacked <- complete(multiimp, 'long') # stack all k=100 imputed datasets
stacked %>% is.na() %>% colSums() # check if all missing values have been imputed

for(i in 1:k){
  imp <- complete(multiimp, i) # extract the ith imputed dataset
  imp$age.cat <- ifelse(imp$RIDAGEYR >= 20, "adult", "adolescent")
  
  imp_log <- full %>% select(SDMVPSU, SDMVSTRA, WTSA2YR) %>% # add back the design variables to the imputed data
    cbind(imp %>% mutate(logPFOA = log(LBXPFOA), # log transformation of PFAS
                         logPFOS = log(LBXPFOS),
                         logPFHxS = log(LBXPFHS),
                         logMeFOSAA = log(LBXMPAH),
                         logPFDA = log(LBXPFDE), 
                         logPFNA = log(LBXPFNA)) ) # %>% 
  # filter(age.cat == "adolescent")
  
  # create a complex survey design object
  design <- svydesign(id = ~ SDMVPSU,     # cluster id: pseudo-PSU
                      strata = ~ SDMVSTRA,  # pseudo-stratum
                      nest = TRUE,       # relabel cluster ids to enforce nesting within strata
                      weights = ~ WTSA2YR,  # sample weights
                      data = imp_log)
  # fit the model (categorical variables have been transformed into factors with reference groups specified)
  # adjusted models (+ RIDAGEYR*edu)
  mod_PFOA <- svyglm(logPFOA ~ carpet.cat4 + RIDAGEYR + gender + race + edu + country + military + INDFMPIR + tapwater + 
                       eatout.cat + shellfish + fish + Standard.Creatinine + period + pregnant.times + children.breastfed, design = design )
  
  mod_PFOS <- svyglm(logPFOS ~ carpet.cat4 + RIDAGEYR + gender + race + edu + country + military + INDFMPIR + tapwater + 
                       eatout.cat + shellfish + fish + Standard.Creatinine + period + pregnant.times + children.breastfed, design = design )
  
  mod_PFHxS <- svyglm(logPFHxS ~ carpet.cat4 + RIDAGEYR + gender + race + edu + country + military + INDFMPIR + tapwater + 
                        eatout.cat + shellfish + fish + Standard.Creatinine + period + pregnant.times + children.breastfed, design = design )
  
  mod_MeFOSAA <- svyglm(logMeFOSAA ~ carpet.cat4 + RIDAGEYR + gender + race + edu + country + military + INDFMPIR + tapwater + 
                          eatout.cat + shellfish + fish + Standard.Creatinine + period + pregnant.times + children.breastfed, design = design )
  
  mod_PFDA <- svyglm(logPFDA ~ carpet.cat4 + RIDAGEYR + gender + race + edu + country + military + INDFMPIR + tapwater + 
                       eatout.cat + shellfish + fish + Standard.Creatinine + period + pregnant.times + children.breastfed, design = design )
  
  mod_PFNA <- svyglm(logPFNA ~ carpet.cat4 + RIDAGEYR + gender + race + edu + country + military + INDFMPIR + tapwater + 
                       eatout.cat + shellfish + fish + Standard.Creatinine + period + pregnant.times + children.breastfed, design = design )
  
  # crude models
  crude_PFOA <- svyglm(logPFOA ~ carpet.cat4, design = design )
  crude_PFOS <- svyglm(logPFOS ~ carpet.cat4, design = design )
  crude_PFHxS <- svyglm(logPFHxS ~ carpet.cat4, design = design )
  crude_MeFOSAA <- svyglm(logMeFOSAA ~ carpet.cat4, design = design )
  crude_PFDA <- svyglm(logPFDA ~ carpet.cat4, design = design )
  crude_PFNA <- svyglm(logPFNA ~ carpet.cat4, design = design )
  
  
  # statistics for adjusted models
  betas1[i,] <- mod_PFOA$coefficients  # extract the betas for the model of logPFOA
  cov.PFOA[,,i] <- vcov(mod_PFOA) # extract the covariance matrix 
  ses1[i,] <- diag(vcov(mod_PFOA))  # extract the variances
  
  betas2[i,] <- mod_PFOS$coefficients  # PFOS
  cov.PFOS[,,i] <- vcov(mod_PFOS)
  ses2[i,] <- diag(vcov(mod_PFOS))
  
  betas3[i,] <- mod_PFHxS$coefficients  
  cov.PFHxS[,,i] <- vcov(mod_PFHxS)
  ses3[i,] <- diag(vcov(mod_PFHxS))
  
  betas4[i,] <- mod_MeFOSAA$coefficients  
  cov.MeFOSAA[,,i] <- vcov(mod_MeFOSAA)
  ses4[i,] <- diag(vcov(mod_MeFOSAA))
  
  betas5[i,] <- mod_PFDA$coefficients  
  cov.PFDA[,,i] <- vcov(mod_PFDA) 
  ses5[i,] <- diag(vcov(mod_PFDA))
  
  betas6[i,] <- mod_PFNA$coefficients  
  cov.PFNA[,,i] <- vcov(mod_PFNA) 
  ses6[i,] <- diag(vcov(mod_PFNA))
  
  
  # statistics for crude models
  betas1.c[i,] <- crude_PFOA$coefficients  # extract the betas for the model of logPFOA
  ses1.c[i,] <- diag(vcov(crude_PFOA))  # extract the variances
  
  betas2.c[i,] <- crude_PFOS$coefficients  # PFOS
  ses2.c[i,] <- diag(vcov(crude_PFOS))
  
  betas3.c[i,] <- crude_PFHxS$coefficients  
  ses3.c[i,] <- diag(vcov(crude_PFHxS))
  
  betas4.c[i,] <- crude_MeFOSAA$coefficients  
  ses4.c[i,] <- diag(vcov(crude_MeFOSAA))
  
  betas5.c[i,] <- crude_PFDA$coefficients  
  ses5.c[i,] <- diag(vcov(crude_PFDA))
  
  betas6.c[i,] <- crude_PFNA$coefficients  
  ses6.c[i,] <- diag(vcov(crude_PFNA))
}

# combine the betas and ses using Rubin's rules (Rubin 1987, page 76)
betas <- list(betas1, betas2, betas3, betas4, betas5, betas6)
ses <- list(ses1, ses2, ses3, ses4, ses5, ses6)

cov1 <- list()
cov1[[1]] <- apply(cov.PFOA, 1:2, mean)
cov1[[2]] <- apply(cov.PFOS, 1:2, mean)
cov1[[3]] <- apply(cov.PFHxS, 1:2, mean)
cov1[[4]] <- apply(cov.MeFOSAA, 1:2, mean)
cov1[[5]] <- apply(cov.PFDA, 1:2, mean)
cov1[[6]] <- apply(cov.PFNA, 1:2, mean)

betas.c <- list(betas1.c, betas2.c, betas3.c, betas4.c, betas5.c, betas6.c)
ses.c <- list(ses1.c, ses2.c, ses3.c, ses4.c, ses5.c, ses6.c)

beta <- cov <- list()

results <- list()

for(i in 1:6){
  
  # average of the k=100 complete-data estimates
  beta[[i]] <- colMeans(betas[[i]])
  # average of the k=100 complete-data variances
  var1 <- colMeans(ses[[i]])
  mat <- betas[[i]] - matrix(beta[[i]], ncol = p, nrow = k, byrow = T) # estimates matrix minus the average estimates by row
  mat_T <- t(mat) # transpose
  
  # between imputations covariance matrix 
  cov2 <- mat_T %*% mat /(k-1)
  # variance between (among) the k=100 complete-data estimates
  # var2 <- diag(mat_T %*% mat)/(k-1)  
  var2 <- diag(cov2)
  
  # total covariance
  cov[[i]] <- cov1[[i]] + (1 + 1/k)*cov2
  #var <- var1 + (1 + 1/k)*var2 # total variance
  var <- diag(cov[[i]])
  # relative increase in variance due to non-response
  r_m <- (1 + 1/k)*var2/var1
  # degrees of freedom
  nu <- (k-1)*(1+1/r_m)^2
  # with alpha = 0.05, upper 2.5 percentile point of the student t distribution on nu degrees of freedom
  t_nu <- qt(p = 0.025, df = nu, lower.tail = FALSE)
  
  CI_hi <- beta[[i]] + t_nu*sqrt(var)  # higher limits of confidence intervals
  CI_low <- beta[[i]] - t_nu*sqrt(var)  # lower limits of confidence intervals 
  
  # exponentiate the estimates and CIs
  beta_exp <- (exp(beta[[i]])-1)*100 
  CI_hi_exp <- (exp(CI_hi)-1)*100
  CI_low_exp <- (exp(CI_low)-1)*100
  
  results[[i]] <- cbind(beta_exp, CI_low_exp, CI_hi_exp)  # combine the exponentiated estimates and CIs
  row.names(results[[i]]) <- tidy(mod_PFOA)$term  # add variable or category names
  
}

View(round(results[[1]], 0))
View(round(results[[2]], 0))
View(round(results[[3]], 0))
View(round(results[[4]], 0))
View(round(results[[5]], 0))
View(round(results[[6]], 0))




# check the r function of linearHypothesis
linearHypothesis
methods(linearHypothesis)
getAnywhere(linearHypothesis.default)



L <- list()
L[[1]] <- matrix(c(0,1,0,0,rep(0,19),
                   0,0,1,0,0,rep(0,18),
                   0,0,0,1,0,0,rep(0,17)),
                 nrow=3, ncol=23, byrow=TRUE)
L[[2]] <- matrix(c(rep(0,6),1,rep(0,16),
                   rep(0,7),1,rep(0,15),
                   rep(0,8),1,rep(0,14)),
                 nrow=3, ncol=23, byrow=TRUE)
L[[3]] <- matrix(c(rep(0,9),1,rep(0,13),
                   rep(0,10),1,rep(0,12)),
                 nrow=2,ncol=23,byrow=TRUE)
L[[4]] <- matrix(c(rep(0,14),1,rep(0,8),
                   rep(0,15),1,rep(0,7)),
                 nrow=2,ncol=23,byrow=TRUE)

pval <- matrix(NA, nrow=6, ncol=4)

for(i in 1:6){
  for(j in 1:4){
    value.hyp <- L[[j]] %*% beta[[i]]
    vcov.hyp <- L[[j]] %*% cov[[i]] %*% t(L[[j]])
    chi <- t(value.hyp) %*% solve(vcov.hyp) %*% value.hyp
    df <- nrow(L[[j]])
    pval[i,j] <- pchisq(chi, df=df, lower.tail = FALSE) 
  }
}
rownames(pval) <- c("PFOA", "PFOS", "PFHxS", "MeFOSAA", "PFDA", "PFNA")
colnames(pval) <- c("Combination of carpet and smooth surface = Low pile carpet = Medium to high pile carpet = Smooth Surface", 
                    "Black = Other Race = White = Hispanic", 
                    "College graduate or above = Some college = Less than college",
                    "Community Supply = Other = Spring = Well or rain cistern = Don't drink tap water")


L <- list()
L[[1]] <- matrix(c(0,1,0,0,rep(0,18),
                   0,0,1,0,0,rep(0,17),
                   0,0,0,1,0,0,rep(0,16)),
                 nrow=3, ncol=22, byrow=TRUE)
L[[2]] <- matrix(c(rep(0,6),1,rep(0,15),
                   rep(0,7),1,rep(0,14),
                   rep(0,8),1,rep(0,13)),
                 nrow=3, ncol=22, byrow=TRUE)
L[[3]] <- matrix(c(rep(0,13),1,rep(0,8),
                   rep(0,14),1,rep(0,7)),
                 nrow=2,ncol=22,byrow=TRUE)

pval <- matrix(NA, nrow=6, ncol=3)

for(i in 1:6){
  for(j in 1:3){
    value.hyp <- L[[j]] %*% beta[[i]]
    vcov.hyp <- L[[j]] %*% cov[[i]] %*% t(L[[j]])
    chi <- t(value.hyp) %*% solve(vcov.hyp) %*% value.hyp
    df <- nrow(L[[j]])
    pval[i,j] <- pchisq(chi, df=df, lower.tail = FALSE) 
  }
}

rownames(pval) <- c("PFOA", "PFOS", "PFHxS", "MeFOSAA", "PFDA", "PFNA")
colnames(pval) <- c("Combination of carpet and smooth surface = Low pile carpet = Medium to high pile carpet = Smooth Surface", 
                    "Black = Other Race = White = Hispanic", 
                    "Community Supply = Other = Spring = Well or rain cistern = Don't drink tap water")
View(pval)



results.c <- list() # results of crude models

for(i in 1:6){
  
  # average of the k=100 complete-data estimates
  beta.c <- colMeans(betas.c[[i]])
  # average of the k=100 complete-data variances
  var1.c <- colMeans(ses.c[[i]])
  mat.c <- betas.c[[i]] - matrix(beta.c, ncol = 4, nrow = k, byrow = T) # estimates matrix minus the average estimates by row
  mat_T.c <- t(mat.c) # transpose
  # variance between (among) the k=100 complete-data estimates
  var2.c <- diag(mat_T.c %*% mat.c)/(k-1)
  # total variance
  var.c <- var1.c + (1 + 1/k)*var2.c
  # relative increase in variance due to non-response
  r_m.c <- (1 + 1/k)*var2.c/var1.c
  # degrees of freedom
  nu.c <- (k-1)*(1+1/r_m.c)^2
  # with alpha = 0.05, upper 2.5 percentile point of the student t distribution on nu degrees of freedom
  t_nu.c <- qt(p = 0.025, df = nu.c, lower.tail = FALSE)
  
  CI_hi.c <- beta.c + t_nu.c*sqrt(var.c)  # higher limits of confidence intervals
  CI_low.c <- beta.c - t_nu.c*sqrt(var.c)  # lower limits of confidence intervals 
  
  # exponentiate the estimates and CIs
  beta_exp.c <- exp(beta.c)  
  CI_hi_exp.c <- exp(CI_hi.c)
  CI_low_exp.c <- exp(CI_low.c)
  
  results.c[[i]] <- cbind(beta_exp.c, CI_low_exp.c, CI_hi_exp.c)  # combine the exponentiated estimates and CIs
  row.names(results.c[[i]]) <- tidy(crude_PFOA)$term  # add variable or category names
  
}


View(results.c[[1]])

#as.mira(m1)
#pool(as.mira(m1))




#svyglm_multiimp <- with(full_multiimp, svyglm(log(LBXPFOA) ~ RIAGENDR + RIDAGEYR + RIDRETH1 + edu + INDFMPIR +  carpet.cat4 + 
#                                                eatout.cat + DRD340 + DRD360 + military + country))
#svyglm_pooled <- pool(svyglm_multiimp)
#summary(svyglm_pooled, conf.int = TRUE, conf.level = 0.95)


# without imputation
### model fitting
pfas <- full[,c("LBXPFOA", "LBXPFOS", "LBXPFHS", "LBXMPAH", "LBXPFDE", "LBXPFNA")]
n1=n2=n3=n4=n5=n6=n7=n8=n9=ci=ci1=ci2=ci3=ci4=ci5=ci6=ci7=ci8=ci9=list()

for(i in 1:6){
  full$Y <- pfas[,i]
  full$lgY <- log(pfas[,i])
  
  des <- svydesign(id = ~ SDMVPSU,
                   strata = ~ SDMVSTRA,
                   nest = TRUE,
                   weights = ~ WTSA2YR,
                   data = full)
  
  n1[[i]] <- svyby(~Y, by =~gender, design = des, FUN=unwtd.count)
  n2[[i]] <- svyby(~Y, by =~race, design = des, FUN=unwtd.count)
  n3[[i]] <- svyby(~Y, by =~edu, design = des, FUN=unwtd.count)
  n4[[i]] <- svyby(~Y, by =~eatout.cat, design = des, FUN=unwtd.count)
  n5[[i]] <- svyby(~Y, by =~shellfish, design = des, FUN=unwtd.count)
  n6[[i]] <- svyby(~Y, by =~fish, design = des, FUN=unwtd.count)
  n7[[i]] <- svyby(~Y, by =~carpet.cat4, design = des, FUN=unwtd.count)
  n8[[i]] <- svyby(~Y, by =~country, design = des, FUN=unwtd.count)
  n9[[i]] <- svyby(~Y, by =~military, design = des, FUN=unwtd.count)
  
  ## Weighted number of individuals
  r <- svymean(~lgY, design=des, na.rm = T) 
  r1 <- svyby(~lgY, by = ~gender, design = des, FUN = svymean, na.rm = T)
  r2 <- svyby(~lgY, by = ~race, design = des, FUN = svymean, na.rm = T)
  r3 <- svyby(~lgY, by = ~edu, design = des, FUN = svymean, na.rm = T)
  r4 <- svyby(~lgY, by = ~eatout.cat, design = des, FUN = svymean, na.rm = T)
  r5 <- svyby(~lgY, by = ~shellfish, design = des, FUN = svymean, na.rm = T)
  r6 <- svyby(~lgY, by = ~fish, design = des, FUN = svymean, na.rm = T)
  r7 <- svyby(~lgY, by = ~carpet.cat4, design = des, FUN = svymean, na.rm = T)
  r8 <- svyby(~lgY, by = ~country, design = des, FUN = svymean, na.rm = T)
  r9 <- svyby(~lgY, by = ~military, design = des, FUN = svymean, na.rm = T)
  
  d <- data.frame(GM=round(exp(r[[1]]), 3))
  d1 <- data.frame(demo.cat = r1[, 1], GM = round(exp(r1[, 2]), 3))
  d2 <- data.frame(demo.cat = r2[, 1], GM = round(exp(r2[, 2]), 3))
  d3 <- data.frame(demo.cat = r3[, 1], GM = round(exp(r3[, 2]), 3))
  d4 <- data.frame(demo.cat = r4[, 1], GM = round(exp(r4[, 2]), 3))
  d5 <- data.frame(demo.cat = r5[, 1], GM = round(exp(r5[, 2]), 3))
  d6 <- data.frame(demo.cat = r6[, 1], GM = round(exp(r6[, 2]), 3))
  d7 <- data.frame(demo.cat = r7[, 1], GM = round(exp(r7[, 2]), 3))
  d8 <- data.frame(demo.cat = r8[, 1], GM = round(exp(r8[, 2]), 3))
  d9 <- data.frame(demo.cat = r9[, 1], GM = round(exp(r9[, 2]), 3))
  
  # 95% CI
  ci[[i]] <- cbind(d, round(exp(confint(r,df=degf(des))), 3))
  ci1[[i]] <- cbind(d1, round(exp(confint(r1,df=degf(des))), 3))
  ci2[[i]] <- cbind(d2, round(exp(confint(r2,df=degf(des))), 3))
  ci3[[i]] <- cbind(d3, round(exp(confint(r3,df=degf(des))), 3))
  ci4[[i]] <- cbind(d4, round(exp(confint(r4,df=degf(des))), 3))
  ci5[[i]] <- cbind(d5, round(exp(confint(r5,df=degf(des))), 3))
  ci6[[i]] <- cbind(d6, round(exp(confint(r6,df=degf(des))), 3))
  ci7[[i]] <- cbind(d7, round(exp(confint(r7,df=degf(des))), 3))
  ci8[[i]] <- cbind(d8, round(exp(confint(r8,df=degf(des))), 3))
  ci9[[i]] <- cbind(d9, round(exp(confint(r9,df=degf(des))), 3))
  
  # svyglm[[i]] <- tidy(svyglm( lgY ~ gender + RIDAGEYR + race + edu + INDFMPIR + carpet.cat4 + eatout.cat + shellfish + fish + tapwater + 
  #                              military + country + period + pregnant.times + children.breastfed, design = des ), 
  #                    exponentiate = TRUE, conf.int = T)
  svyglm[[i]] <- svyglm( lgY ~ gender + RIDAGEYR + race + edu + INDFMPIR + carpet.cat4 + eatout.cat + shellfish + fish + tapwater + 
                           military + country + period + pregnant.times + children.breastfed + Standard.Creatinine, design = des )
}
# + military + country
names(n1) <- colnames(pfas)
n1

names(n2) <- colnames(pfas)
n2

names(n3) <- colnames(pfas)
n3

names(n4) <- colnames(pfas)
n4

names(n5) <- colnames(pfas)
n5

names(n6) <- colnames(pfas)
n6

names(n7) <- colnames(pfas)
n7

names(n8) <- colnames(pfas)
n8

names(n9) <- colnames(pfas)
n9

names(ci) <- colnames(pfas)
ci

names(ci1) <- colnames(pfas)
ci1

names(ci2) <- colnames(pfas)
ci2

names(ci3) <- colnames(pfas)
ci3

names(ci4) <- colnames(pfas)
ci4

names(ci5) <- colnames(pfas)
ci5

names(ci6) <- colnames(pfas)
ci6

names(ci7) <- colnames(pfas)
ci7

names(ci8) <- colnames(pfas)
ci8

names(ci9) <- colnames(pfas)
ci9

names(svyglm) <- colnames(pfas)
tidy(svyglm$LBXPFOA, conf.int = TRUE, exponentiate = TRUE) %>% 
  mutate(estimate = (estimate - 1)*100, 
         conf.low = (conf.low - 1)*100, 
         conf.high = (conf.high - 1)*100 ) %>% View()

tidy(svyglm$LBXPFOS, conf.int = TRUE, exponentiate = TRUE) %>% 
  mutate(estimate = (estimate - 1)*100, 
         conf.low = (conf.low - 1)*100, 
         conf.high = (conf.high - 1)*100 ) %>% View()

tidy(svyglm$LBXPFHS, conf.int = TRUE, exponentiate = TRUE) %>% 
  mutate(estimate = (estimate - 1)*100, 
         conf.low = (conf.low - 1)*100, 
         conf.high = (conf.high - 1)*100 ) %>% View()

tidy(svyglm$LBXMPAH, conf.int = TRUE, exponentiate = TRUE) %>% 
  mutate(estimate = (estimate - 1)*100, 
         conf.low = (conf.low - 1)*100, 
         conf.high = (conf.high - 1)*100 ) %>% View()

tidy(svyglm$LBXPFDE, conf.int = TRUE, exponentiate = TRUE) %>% 
  mutate(estimate = (estimate - 1)*100, 
         conf.low = (conf.low - 1)*100, 
         conf.high = (conf.high - 1)*100 ) %>% View()

tidy(svyglm$LBXPFNA, conf.int = TRUE, exponentiate = TRUE) %>% 
  mutate(estimate = (estimate - 1)*100, 
         conf.low = (conf.low - 1)*100, 
         conf.high = (conf.high - 1)*100 ) %>% View()


b <- s <- list()
for(i in 1:6){
  b[[i]] <- svyglm[[i]]$coefficients %>% as.matrix()
  s[[i]] <- vcov(svyglm[[i]]) %>% as.matrix()
}


L1 <- list()

L1[[1]] <- matrix(c(0,0,0,1,rep(0,19),
                    0,0,0,0,1,rep(0,18),
                    0,0,0,0,0,1,rep(0,17)),
                  nrow=3, ncol=23, byrow=TRUE)
L1[[2]] <- matrix(c(rep(0,6),1,rep(0,16),
                    rep(0,7),1,rep(0,15)),
                  nrow=2, ncol=23, byrow=TRUE)
L1[[3]] <- matrix(c(rep(0,9),1,rep(0,13),
                    rep(0,10),1,rep(0,12),
                    rep(0,11),1,rep(0,11)),
                  nrow=3,ncol=23,byrow=TRUE)
L1[[4]] <- matrix(c(rep(0,15),1,rep(0,7),
                    rep(0,16),1,rep(0,6)),
                  nrow=2,ncol=23,byrow=TRUE)

pval.1 <- matrix(NA, nrow=6, ncol=4)

for(i in 1:6){
  for(j in 1:4){
    value.hyp <- L1[[j]] %*% b[[i]]
    vcov.hyp <- L1[[j]] %*% s[[i]] %*% t(L1[[j]])
    chi <- t(value.hyp) %*% solve(vcov.hyp) %*% value.hyp
    df <- nrow(L1[[j]])
    pval.1[i,j] <- pchisq(chi, df=df, lower.tail = FALSE) 
  }
}

rownames(pval.1) <- c("PFOA", "PFOS", "PFHxS", "MeFOSAA", "PFDA", "PFNA")
colnames(pval.1) <- c("Black = Other Race = White = Hispanic", 
                      "College graduate or above = Some college = Less than college",
                      "Combination of carpet and smooth surface = Low pile carpet = Medium to high pile carpet = Smooth Surface", 
                      "Community Supply = Spring = Well or rain cistern = Don't drink tap water")
View(pval.1)




PFAS_ALDUST_DEMO_2005 %>% filter(!is.na(carpet.cat)) %>% filter(!is.na(edu)) %>% 
  filter(!is.na(DRD340)) %>% filter(!is.na(eatout.cat)) %>% filter(!is.na(DRD360)) %>% View()

PFAS_ALDUST_DEMO_2005 %>% filter(!is.na(carpet.cat)) %>% filter(!is.na(edu)) %>% 
  filter(!is.na(DRD340)) %>% filter(!is.na(eatout.cat)) %>% filter(!is.na(DRD360)) %>% 
  filter(DR1DRSTZ == 2) %>% View()

PFAS_ALDUST_DEMO_2005 %>% filter(!is.na(gender) & !is.na(RIDAGEYR) & !is.na(race)) %>%
  filter(!is.na(edu) & !is.na(INDFMPIR) & !is.na(carpet.cat4) & !is.na(eatout.cat)) %>%
  filter(!is.na(shellfish) & !is.na(fish)) %>% filter(!is.na(military) & !is.na(country)) %>% 
  filter(!is.na(tapwater) & !is.na(period) & !is.na(pregnant.times) & !is.na(children.breastfed) ) %>% 
  nrow() #1108

full %>% filter(!is.na(gender) & !is.na(RIDAGEYR) & !is.na(race)) %>%
  filter(!is.na(edu) & !is.na(INDFMPIR) & !is.na(carpet.cat4) & !is.na(eatout.cat)) %>%
  filter(!is.na(shellfish) & !is.na(fish)) %>% filter(!is.na(military) & !is.na(country)) %>% 
  filter(!is.na(tapwater) & !is.na(period) & !is.na(pregnant.times) & !is.na(children.breastfed) & !is.na(LBXSCR) ) %>% 
  nrow() # 1044


PFAS_ALDUST_DEMO_2005.1 <- PFAS_ALDUST_DEMO_2005 %>% filter(DR1DRSTZ == 1)
PFAS_ALDUST_DEMO_2005.1 %>% filter(!is.na(RIAGENDR) & !is.na(RIDAGEYR) & !is.na(RIDRETH1)) %>%
  filter(!is.na(edu) & !is.na(INDFMPIR) & !is.na(carpet.cat4) & !is.na(eatout.cat)) %>%
  filter(!is.na(DRD340) & !is.na(DRD360)) %>% nrow() #1413


n = ci1 = ci2 = svyglm = list()

for(i in 1:6){
  PFAS_ALDUST_DEMO_2005$Y <- pfas[,i]
  PFAS_ALDUST_DEMO_2005$lgY <- log(pfas[,i])
  
  des <- svydesign(id = ~ SDMVPSU,
                   strata = ~ SDMVSTRA,
                   nest = TRUE,
                   weights = ~ WTSA2YR,
                   data = PFAS_ALDUST_DEMO_2005)
  
  n[[i]] <- svyby(~Y, by=~carpet.cat4, design = des, FUN=unwtd.count)
  
  ## Weighted number of individuals
  r1 <- svymean(~lgY, design=des, na.rm = T) 
  r2 <- svyby(~lgY, by = ~carpet.cat4, design = des, FUN = svymean, na.rm = T)
  
  d1 <- data.frame(GM=round(exp(r1[[1]]), 3))
  d2 <- data.frame(demo.cat = r2[, 1], GM = round(exp(r2[, 2]), 3))
  
  # 95% CI
  ci1[[i]] <- cbind(d1, round(exp(confint(r1,df=degf(des))), 3))
  ci2[[i]] <- cbind(d2, round(exp(confint(r2,df=degf(des))), 3))
  
  svyglm[[i]] <- tidy(svyglm( lgY ~ as.factor(carpet.cat4), design = des ), 
                      exponentiate = TRUE, conf.int = T)
  
}

names(n) <- colnames(pfas)
n

names(ci1) <- colnames(pfas)
ci1

names(ci2) <- colnames(pfas)
ci2

names(svyglm) <- colnames(pfas)
svyglm




des <- svydesign(id = ~ SDMVPSU,
                 strata = ~ SDMVSTRA,
                 nest = TRUE,
                 weights = ~ WTSA2YR,
                 data = PFAS_ALDUST_DEMO_2005)
svychisq(~ RIAGENDR + military, design = des, statistics = "Chisq")
svychisq(~ carpet.cat + country, design = des, statistics = "Chisq")
svychisq(~ edu + country, design = des, statistics = "Chisq")
svychisq(~ eatout.cat + country, design = des, statistics = "Chisq")
svychisq(~ military + country, design = des, statistics = "Chisq")
svychisq(~ RIDRETH1 + country, design = des, statistics = "Chisq")

library(reshape2)
corr <- round(cor(PFAS_ALDUST_DEMO_2005.2 %>% 
                    select(LBXPFOA, LBXPFOS, LBXPFHS, LBXMPAH, LBXPFDE, LBXPFNA) %>% 
                    filter(!is.na(LBXPFOA)) %>% 
                    rename(PFOA = LBXPFOA, PFOS = LBXPFOS, PFHxS = LBXPFHS, MeFOSAA = LBXMPAH, PFDA = LBXPFDE, PFNA = LBXPFNA)), 2)
melted.corr <- melt(corr)

dat <- PFAS_ALDUST_DEMO_2005 %>% filter(!is.na(LBXPFOA)) %>%
  select(SDMVPSU, SDMVSTRA, WTSA2YR, LBXPFOA, LBXPFOS, LBXPFHS, LBXMPAH, LBXPFDE, LBXPFNA) %>%
  mutate(PFOA = dense_rank(LBXPFOA), PFOS = dense_rank(LBXPFOS), PFHxS = dense_rank(LBXPFHS), 
         MeFOSAA = dense_rank(LBXMPAH), PFDA = dense_rank(LBXPFDE), PFNA = dense_rank(LBXPFNA))

library(survey)
des.1 <- svydesign(id = ~ SDMVPSU,
                   strata = ~ SDMVSTRA,
                   nest = TRUE,
                   weights = ~ WTSA2YR,
                   data = dat)
library(jtools)
cor <- svycor(~ PFOA + PFOS + PFHxS + MeFOSAA + PFDA + PFNA, design = des.1)
melted.cor <- melt(cor$cors)

### Spearman Correlation Heatmap

tiff("Heatmap.tiff", width = 5, height = 5, units = 'in', res = 600)

library(ggplot2)
ggplot(data = melted.cor, aes(Var1, Var2, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Spearman\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, 
                                   hjust = 1)) + xlab("") + ylab("") +
  coord_fixed()
dev.off()

library(ggcorrplot)

tiff("Heatmap.1.tiff", width = 5, height = 5, units = 'in', res = 600)
ggcorrplot(cor[["cors"]], hc.order = TRUE, outline.col = "white", lab = TRUE, legend.title = "Spearman\nCorrelation")
dev.off()

### forest plot
workdir <- " "
Carpetingfile <- file.path(workdir,"Carpeting.csv")
Carpeting <- read.csv(Carpetingfile, stringsAsFactors=FALSE)

Carpetingtext <- cbind(c("PFAS", "\n", Carpeting$X), 
                       c("Model", "\n", Carpeting$X.1))

library(forestplot)
png(file.path(workdir,"PFASplot.png"),width=240, height=135)
tiff(file.path(workdir, "Figures/PFASplot.tiff"), width = 12, height = 8, units = 'in', res = 600)
forestplot(labeltext=Carpetingtext, graph.pos=3, 
           mean=c(NA, NA, Carpeting$Point.Estimate), 
           lower=c(NA, NA, Carpeting$Low), upper=c(NA, NA, Carpeting$High),
           title="Confidence Interval (ng/ml)",
           hrzl_lines=list("3" = gpar(lwd=1, col="#99999922")),
           xticks=c(-1, 0, 8),
           txt_gp=fpTxtGp(label=gpar(cex=1.25),
                          ticks=gpar(cex=1.1),
                          xlab=gpar(cex = 1.2),
                          title=gpar(cex = 1.2)),
           col=fpColors(box="#1c61b6", lines="#1c61b6", zero = "gray50"),
           zero=0, cex=0.9, lineheight = "auto", boxsize=0.3, colgap=unit(10,"mm"),
           lwd.ci=2, ci.vertices=TRUE, ci.vertices.height = 0.2)

dev.off()


PFAS_ALDUST_DEMO_2005.3 <- PFAS_ALDUST_DEMO_2005.2 %>% left_join(BMX_D.1, by = "SEQN")



### Boxplot 
library(ggplot2)
ggplot(data = full, aes(x = carpet.cat4, y = LBXPFOA)) + geom_boxplot(aes(fill = carpet.cat4)) +
  theme(axis.text.x = element_text(angle = 45))
ggplot(data = full, aes(x = carpet.cat4, y = LBXPFOS)) + geom_boxplot(aes(fill = carpet.cat4)) +
  theme(axis.text.x = element_text(angle = 45))
ggplot(data = full, aes(x = carpet.cat4, y = LBXPFHS)) + geom_boxplot(aes(fill = carpet.cat4)) +
  theme(axis.text.x = element_text(angle = 45))


### plot
data <- read_excel("data.xlsx")

ggplot(data = data, aes(x = Type, y = GeometricMean, ymin = LowerLimit, ymax = UpperLimit )) +
  geom_pointrange(aes(col = Type)) +
  xlab('Type of Flooring') + ylab("Geometric Mean (95% Confidence Interval)") +
  geom_errorbar(aes(ymin = LowerLimit, ymax = UpperLimit, col = Type),width = 0.5, cex = 1) + 
  facet_wrap(~ PFAS, strip.position = "left", nrow = 9, scales = "free_y") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(face = "bold"),
        axis.title = element_text(size = 12, face = "bold"),
        strip.text.y = element_text(hjust = 0, vjust = 1, angle = 180, face = "bold"))+
  coord_flip() +
  scale_y_log10() 

data$PFAS <- factor(data$PFAS, levels = c('PFOA', 'PFOS', 'PFHxS', 'MeFOSAA', 'PFDA', 'PFNA'))

ggplot(data = data, aes(x = Type, y = GeometricMean, ymin = LowerLimit, ymax = UpperLimit )) +
  geom_pointrange(aes(col = Type)) +
  xlab('Type of Flooring') + ylab("Geometric Mean (95% Confidence Interval)") +
  geom_errorbar(aes(ymin = LowerLimit, ymax = UpperLimit, col = Type), width = 0.5, cex = 1) + 
  facet_wrap(~ PFAS, strip.position = "top", nrow = 1, scales = "free_x") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 12, face = "bold")) +
  scale_colour_discrete("Type of Flooring") + 
  scale_y_log10() 


figure1 <- 
  ggplot(data = data %>% filter(PFAS == "PFOA"), aes(x = Type, y = GeometricMean, ymin = LowerLimit, ymax = UpperLimit)) +
  geom_pointrange(aes(col = Type)) + 
  xlab("") + ylab("") + 
  geom_errorbar(aes(ymin = LowerLimit, ymax = UpperLimit, col = Type), width = 0.5, cex = 1) +
  facet_wrap(~ PFAS, strip.position = "top", nrow = 1, scales = "free_x") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 12, face = "bold")) +
  scale_colour_discrete("Type of Flooring") 


figure2 <-
  ggplot(data = data %>% filter(PFAS == "PFOS"), aes(x = Type, y = GeometricMean, ymin = LowerLimit, ymax = UpperLimit)) +
  geom_pointrange(aes(col = Type)) + 
  xlab("") + ylab("") +
  geom_errorbar(aes(ymin = LowerLimit, ymax = UpperLimit, col = Type), width = 0.5, cex = 1) +
  facet_wrap(~ PFAS, strip.position = "top", nrow = 1, scales = "free_x") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 12, face = "bold")) +
  scale_colour_discrete("Type of Flooring") 

figure3 <-
  ggplot(data = data %>% filter(PFAS == "PFHxS"), aes(x = Type, y = GeometricMean, ymin = LowerLimit, ymax = UpperLimit)) +
  geom_pointrange(aes(col = Type)) + 
  xlab("") + ylab("") +
  geom_errorbar(aes(ymin = LowerLimit, ymax = UpperLimit, col = Type), width = 0.5, cex = 1) +
  facet_wrap(~ PFAS, strip.position = "top", nrow = 1, scales = "free_x") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 12, face = "bold")) +
  scale_colour_discrete("Type of Flooring") 

figure4 <-
  ggplot(data = data %>% filter(PFAS == "MeFOSAA"), aes(x = Type, y = GeometricMean, ymin = LowerLimit, ymax = UpperLimit)) +
  geom_pointrange(aes(col = Type)) + 
  xlab("") + ylab("") + 
  geom_errorbar(aes(ymin = LowerLimit, ymax = UpperLimit, col = Type), width = 0.5, cex = 1) +
  facet_wrap(~ PFAS, strip.position = "top", nrow = 1, scales = "free_x") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 12, face = "bold")) +
  scale_colour_discrete("Type of Flooring") 

figure5 <-
  ggplot(data = data %>% filter(PFAS == "PFDA"), aes(x = Type, y = GeometricMean, ymin = LowerLimit, ymax = UpperLimit)) +
  geom_pointrange(aes(col = Type)) + 
  xlab('') + ylab("") +
  geom_errorbar(aes(ymin = LowerLimit, ymax = UpperLimit, col = Type), width = 0.5, cex = 1) +
  facet_wrap(~ PFAS, strip.position = "top", nrow = 1, scales = "free_x") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 12, face = "bold")) +
  scale_colour_discrete("Type of Flooring") 

figure6 <-
  ggplot(data = data %>% filter(PFAS == "PFNA"), aes(x = Type, y = GeometricMean, ymin = LowerLimit, ymax = UpperLimit)) +
  geom_pointrange(aes(col = Type)) + 
  xlab("") + ylab("") +
  geom_errorbar(aes(ymin = LowerLimit, ymax = UpperLimit, col = Type), width = 0.5, cex = 1) +
  facet_wrap(~ PFAS, strip.position = "top", nrow = 1, scales = "free_x") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 12, face = "bold")) +
  scale_colour_discrete("Type of Flooring") 


workdir <- ""
tiff(file.path(workdir, "PFASGeometricMeans.tiff"), width = 16, height = 8, units = 'in', res = 600)
figure <- ggarrange(figure1, figure2, figure3, figure4, figure5, figure6,
                    labels = c("A", "B", "C", "D", "E", "F"),
                    ncol = 3, nrow = 2)
# Annotate the figure by adding a common labels
figure <- 
  annotate_figure(figure,
                  #top = text_grob("Visualizing len", color = "red", face = "bold", size = 14),
                  bottom = text_grob("Type of Flooring", color = "black", size = 12, face = "bold"),
                  left = text_grob("Geometric Mean (95% Confidence Interval)", color = "black", rot = 90, size = 12, face = "bold"))
figure
dev.off()


library(readxl)
data1 <- read_excel("data1.xlsx")

data1$PFAS <- factor(data1$PFAS, levels = c('PFOA', 'PFOS', 'PFHxS', 'MeFOSAA', 'PFDA', 'PFNA'))

tiff(file.path(workdir, "MissingData5.tiff"), width = 12, height = 7, units = 'in', res = 600)

ggplot(data = data1, aes(x = Type, y = `% difference in serum PFAS concentration`, ymin = LowerLimit, ymax = UpperLimit )) +
  geom_pointrange(aes(col = Type)) +
  geom_hline(aes(fill = Type), yintercept = 1, linetype = 2) +
  xlab('Type of Flooring') + ylab("% difference in serum PFAS concentration (95% Confidence Interval)") +
  geom_errorbar(aes(ymin = LowerLimit, ymax = UpperLimit, col = Type), width = 0.5, cex = 1) + 
  facet_wrap(~ PFAS, strip.position = "top", nrow = 1, scales = "free_x") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 12, face = "bold")) +
  scale_colour_discrete("Type of Flooring")

dev.off()


### complete case analysis results
data2 <- read_excel("data2.xlsx")

data2$PFAS <- factor(data2$PFAS, levels = c('PFOA', 'PFOS', 'PFHxS', 'MeFOSAA', 'PFDA', 'PFNA'))

tiff(file.path(workdir, "MissingData6.tiff"), width = 12, height = 7, units = 'in', res = 600)

ggplot(data = data2, aes(x = Type, y = `% difference in serum PFAS concentration`, ymin = LowerLimit, ymax = UpperLimit )) +
  geom_pointrange(aes(col = Type)) +
  geom_hline(aes(fill = Type), yintercept = 1, linetype = 2) +
  xlab('Type of Flooring') + ylab("% difference in serum PFAS concentration (95% Confidence Interval)") +
  geom_errorbar(aes(ymin = LowerLimit, ymax = UpperLimit, col = Type), width = 0.5, cex = 1) + 
  facet_wrap(~ PFAS, strip.position = "top", nrow = 1, scales = "free_x") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 12, face = "bold")) +
  scale_colour_discrete("Type of Flooring")

dev.off()


### Adults
data3 <- read_excel("data3.xlsx")

data3$PFAS <- factor(data3$PFAS, levels = c('PFOA', 'PFOS', 'PFHxS', 'MeFOSAA', 'PFDA', 'PFNA'))

tiff(file.path(workdir, "Adults.tiff"), width = 12, height = 4, units = 'in', res = 600)

ggplot(data = data3, aes(x = Type, y = `% difference in serum PFAS concentration`, ymin = LowerLimit, ymax = UpperLimit )) +
  geom_pointrange(aes(col = Type)) +
  geom_hline(aes(fill = Type), yintercept = 1, linetype = 2) +
  xlab('Type of Flooring') + ylab("% difference in serum PFAS concentration (95% Confidence Interval)") +
  geom_errorbar(aes(ymin = LowerLimit, ymax = UpperLimit, col = Type), width = 0.5, cex = 1) + 
  facet_wrap(~ PFAS, strip.position = "top", nrow = 1, scales = "free_x") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 12, face = "bold")) +
  scale_colour_discrete("Type of Flooring")

dev.off()


### Adolescents
data4 <- read_excel("data4.xlsx")

data4$PFAS <- factor(data4$PFAS, levels = c('PFOA', 'PFOS', 'PFHxS', 'MeFOSAA', 'PFDA', 'PFNA'))

tiff(file.path(workdir, "Adolescents.tiff"), width = 12, height = 4, units = 'in', res = 600)

ggplot(data = data4, aes(x = Type, y = `% difference in serum PFAS concentration`, ymin = LowerLimit, ymax = UpperLimit )) +
  geom_pointrange(aes(col = Type)) +
  geom_hline(aes(fill = Type), yintercept = 1, linetype = 2) +
  xlab('Type of Flooring') + ylab("% difference in serum PFAS concentration (95% Confidence Interval)") +
  geom_errorbar(aes(ymin = LowerLimit, ymax = UpperLimit, col = Type), width = 0.5, cex = 1) + 
  facet_wrap(~ PFAS, strip.position = "top", nrow = 1, scales = "free_x") +
  theme(plot.title = element_text(size = 16, face = "bold"),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 12, face = "bold")) +
  scale_colour_discrete("Type of Flooring")

dev.off()


### forest plot for all variables
library(forestplot)
data6 <- read_excel("data6_PFOA.xlsx")

data6$Variable <- c("Type of floor covering", "   Smooth surface (ref)", "   Low pile carpet",  "   Mixed surface", 
                    "   Medium to high pile carpet", NA, "Age", NA, "Gender", "   Female (ref)", "   Male", NA, 
                    "Race/ethnicity", "   Hispanic (ref)", "   Non-Hispanic White", "   Non-Hispanic Black", "   Other Race",  
                    NA, "Education", "   Less than college (ref)", "   Some college", "   College graduate or above", NA, 
                    "Country of birth", "   Foreign (ref)", "   Born in the US", NA, "Served in the military", 
                    "   No (ref)", "   Yes", NA, "Family PIR", NA, "Tap water source", "   Don't drink tap water (ref)",
                    "   Community supplied water", "   Other water source", NA, "Eatout every week", 
                    "   No (ref)", "   Yes", NA, "Had shellfish in past 30 days", "   No (ref)", "   Yes", NA, 
                    "Had fish in past 30 days", "   No (ref)", "   Yes", NA, "Standard creatinine", NA, 
                    "Had at least one period in the past one year", "    No (ref)", "    Yes", NA, 
                    "Number of pregnancies", NA, "Number of children breastfed")                

plottext <- c("\n", "\n", data6$Variable)

tiff(file.path("PFOAplot1.tiff"), width = 12, height = 18, units = 'in', res = 600)
forestplot(labeltext = plottext, graph.pos = 2, 
           mean=c(NA, NA, data6$beta), 
           lower=c(NA, NA, data6$CI_low_exp), upper=c(NA, NA, data6$CI_hi_exp),
           title="Confidence Interval (ng/ml)",
           hrzl_lines=list("3" = gpar(lwd=1, col="#99999922")),
           xticks=c(-45, 0, 45),
           txt_gp=fpTxtGp(label=gpar(cex=1.25),
                          ticks=gpar(cex=1.1),
                          xlab=gpar(cex = 1.2),
                          title=gpar(cex = 1.2)),
           col=fpColors(box="#1c61b6", lines="#1c61b6", zero = "gray50"),
           zero=0, cex=0.9, lineheight = "auto", boxsize=0.3, colgap=unit(10,"mm"),
           lwd.ci=2, ci.vertices=TRUE, ci.vertices.height = 0.2)

dev.off()
