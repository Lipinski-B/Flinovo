### Kaplan Meyer (categorical variables)
### COX (continuous variables)

library(survival)
library(survminer)
library(ggplot2)
library(ggpubr)
library(ggfortify)
library(ranger)
library(ggpubr)
library(grid)
library(readxl)
library(ggpubr)
library(forestplot)
library(forestmodel)

# Cohort censored at 10 years = 3650 days

# FLIPI score 0 ou 1 = 1_Low
# FLIPI score 2 = 2_Inter
# FLIPI score >=3 3_High

# PRIMA score 1_Low = 0
# PRIMA score Int = 1 = 2_High

# Load dataframe df2_Fig1G_S1H.xlsx

df2[df2 == "NA"] <- NA

### Kaplan Meier analysis

km_fit <- survfit(Surv(PFS, Relapse) ~ 1, data=df2)
plot <- ggsurvplot(km_fit, pval = TRUE, conf.int = FALSE, palette = "blue", risk.table = TRUE, 
                   xlab = "Time (months)", ylab = "Progression free survival probability", ggtheme = theme_classic(), legend = "none", tables.theme = theme_cleantable(), risk.table.y.text = FALSE,
                   fontsize = 4)
plot$table <- plot$table + theme(plot.title = element_text(hjust = -0.05)) + theme(plot.title = element_text(size=10))
plot


### COX

cox_Red <- coxph(Surv(PFS, Relapse) ~ Volume_reduction, data = df2)

cox_PRIMA <- coxph(Surv(PFS, Relapse) ~ PRIMA_PI, data = df2)

cox_FLIPI <- coxph(Surv(PFS, Relapse) ~ FLIPI, data = df2)

cx1 <- coxph(Surv(PFS, Relapse) ~ Volume_reduction + FLIPI, data = as.data.frame(df2)) 

forest_model(cx1, covariates = c("Volume_reduction", "FLIPI"), format_options = forest_model_format_options(text_size = 4, point_size = 2), recalculate_width = TRUE)
