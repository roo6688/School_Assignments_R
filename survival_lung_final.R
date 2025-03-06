library(survival)
library(survminer)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(GGally)

### Chargement du jeu de données "lung" ###
lung
lung <- data.frame(cancer, package="survival")

#### Analyse exploratoire des données ####
summary(lung)

# Correlations entre chacune des variables 
ggpairs(lung, columns = (2:9), upper = list(continuous = wrap("cor", size = 4)))

# Correlations entre le temps de survie et le score ECOG (time~ph.ecog)
lung %>% 
  ggplot(aes(time, ph.ecog)) + 
  geom_point() + 
  geom_smooth(method = "lm")
  
lm_ph.ecog <- lm(time~ph.ecog, lung)
summary(lm_ph.ecog)
boxplot(time~ph.ecog, data = lung,
        xlab = "Score ECOG",
        ylab = "Jours")

##### Courbes de survie Kaplan-Meier ####
# Overall fit
fit_overall = survfit(Surv(time, status) ~ 1, data=lung) #BASE R
plot(fit_overall,
     xlab = "Jours",
     ylab = "Probabilité de survie")
print(fit_overall)

ggsurvplot(fit_overall,  # SURVMINER
           xlab="Jours", 
           ylab="Probabilité de survie",
           legend = c(0.8,0.8),
           legend.title = "",
           legend.labs = "Tout",
           risk.table="abs_pct",
           risk.table.title = "Nombre et % de patients à risque",
           risk.table.labs = "",
           risk.table.height = 0.18,
           conf.int=TRUE,
           surv.median.line="hv",
           cumcensor = TRUE,
           cumcensor.title = "Nombre cumulatif de censures",
           cumcensor.labs = "",
           cumcensor.height = 0.18,
           censor.size = 7)

# ph.ecog
km_fit1 <- survfit(Surv(time, status) ~ ph.ecog, data = lung)

plot(km_fit1, col = c("blue", "red", "green", "purple"), # BASE R
     mark.time=TRUE, pch=8, 
     xlab="Jours", ylab="Probabilité de survie",
     lwd = 1.25)
legend("topright", title = "Score ECOG évalué par médecin", 
       legend = c("0 = asymptomatique",
                  "1 = symptomatique mais ambulatoire",
                  "2 = au lit < 50% de la journée",
                  "3 - au lit > 50 % de la journée mais pas alité"),
       col = c("blue", "red", "green", "purple"), cex = 0.55, lwd = 1.25)

summary(lung$ph.ecog)
count(lung,ph.ecog)

# ph.ecog (sans ECOG=3, plot survminer)
lung2 <- lung %>% 
  filter(ph.ecog != 3) # Pour enlever le outlier

km_fit2 <- survfit(Surv(time, status) ~ ph.ecog, data = lung2) #SURVMINER
ggsurvplot(km_fit2, 
           xlab="Jours", 
           ylab="Probabilité de survie",
           main = "Survie par score de performance ECOG",
           legend.title = "Score ECOG",
           legend.labs = c("0","1","2"),
           risk.table="abs_pct",
           risk.table.title = "Nombre et % de patients à risque",
           conf.int=TRUE,
           surv.median.line="hv",
           cumcensor = TRUE,
           cumcensor.title = "Nombre cumulatif de censures",
           censor.size = 7)
count(lung2,ph.ecog)
km_fit2

#### Cox Proportional Hazards Model ####
# ph.ecog (univariate)
cox_fit1 <- coxph(Surv(time, status) ~ ph.ecog, data = lung2)
summary(cox_fit1)

#ggforest(cox_fit1) ?

# ph.ecog + wt.loss (multivariate)
cox_fit2 <- coxph(Surv(time, status) ~ ph.ecog + wt.loss, data = lung2)
summary(cox_fit2)

# Évaluation de performance (ph.ecog + ph.karno + pat.karno, multivariate) -- EXTRA
cox_fit3 <- coxph(Surv(time, status) ~ ph.ecog + ph.karno + pat.karno, data = lung2)
summary(cox_fit3)

#### LOG RANK TEST #####
# ph.ecog
log_rank <- survdiff(Surv(time, status) ~ ph.ecog, lung2)
log_rank

