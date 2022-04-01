# Estatistika -------------------------------------------------------------

library(ggplot2)
library(plyr)

x <- openxlsx::read.xlsx("G:/Mi unidad/NASIR/MAIN_FOLDER_NASIR/MAIN_TABLE_v014_newPFS_TTP.xlsx", 2)

library(vtable)
st(x)

library(gtsummary)
# make dataset with a few variables to summarize
trial2 <- x %>% select(BENEFIT2_39, Age_34, T_D_NAG_Strong500RNA1_227, Idibell_D_SigX16_272, T_IF_Tregs_PD1_pr_Tumor_168)

# summarize the data with our package
table1 <- tbl_summary(trial2, by = BENEFIT2_39) %>% add_p() %>% add_overall() %>% bold_labels()  %>% add_n()

mydata <- summarytools::descr(x)
# View(mydata)


lehen_var <- grep("OS_days_68", names(x))
new_fact <- c("Age_34", names(x)[lehen_var:ncol(x)])
constraste <- c("BENEFIT2_39","RESPUESTA_43","CONTROL_44","Progresion_45",
                "PVI_29","TACE_30","Sorafenib_31","Etiology_32","AgeMay65_33","AFPMzy400_35","MAAuptake_36",
                "PFS_49","Muerte_51","Resection_53","Nivolumabpost_54","TKIpost_55", "CRPR_SDPD_bin_61b",
                "PatronPD_en2_12_48b")

for(c in constraste)
{
  print(c)
  x[,c] <- as.factor(x[,c])
}


# Contraste simple --------------------------------------------------------



normal_df <- NULL

for(h in new_fact)
{
  # h <- new_fact[1]
  
  for(co in constraste)
  {
    # co <- constraste[4]
    print(paste0("Contrast ", co, " - Variable ", h))
    temp_df <- x[,c(co, h)]
    temp2 <- temp_df[complete.cases(temp_df),]
    
    if(nrow(temp2) > 1 & table(temp2[,co])[1] > 1 & table(temp2[,co])[2] > 1 & var(temp2[,h]) != 0)
    {
      print(paste0("Contrast ", co, " - Variable ", h, "-----------------> SI"))
      temp2[,co] <- factor(temp2[,co])
      # mu <- ddply(temp2, co, summarise, grp.mean=mean(temp2[,h]))
      mu <- aggregate(temp2[,h], by = list(temp2[,co]), mean)
      names(mu) <- c(co, "grp_mean")
      
      # plot
      # p1 <- ggplot(temp2, aes_string(x=h)) +
      #   geom_density(alpha=0.4)
      #
      # p2 <- ggplot(temp2, aes_string(x=h, fill=co)) +
      #   geom_density(alpha=0.4) +geom_vline(data=mu, aes(xintercept=grp_mean, color=mu[,co]),
      #                                       linetype="dashed")
      #
      #
      # jpeg(paste0("Images/Fac_",h,"_contrs_",co,"_general.jpg"))
      # print(p1)
      # dev.off()
      #
      # jpeg(paste0("Images/Fac_",h,"_contrs_",co,"_bygroup.jpg"))
      # print(p2)
      # dev.off()
      
      # normal?
      # Bigger than 0.05 > p-value --> normal
      st <- NULL
      
      if(nrow(temp2) > 3)
      {
        st <- shapiro.test(temp2[,h])
        
        normal <- ifelse(st$p.value > 0.05, "Normal","No-Normal")
        
        if(nrow(temp2) > 4)
        {
          
          if(normal == "Normal")
          {
            ts <- t.test(temp2[,h] ~ temp2[,co])
          } else {
            ts <- wilcox.test(temp2[,h] ~ temp2[,co])
          }
          
          important <- ifelse(ts$p.value > 0.05, "x","Important")
          
          ag1 <- aggregate(temp2[,h], by = list(temp2[,co]), mean)
          
          normal_df <- c(normal_df, h, co, nrow(temp2),  table(temp2[,co])[1],  table(temp2[,co])[2], 
                         round(ag1$x[1],2), round(ag1$x[2],2),
                         st$p.value, normal, ts$p.value, important)
          
        }
      }
    }
    
    
  }
}

df_normal <- data.frame(matrix(normal_df, ncol = 11, byrow = TRUE))
names(df_normal) <- c("Variable","Contrast","Total patients","Group0","Group1",
                      "Mean0","Mean1",
                      "NormalPvalue","Normal","p-value","Important")

head(df_normal)

for(n in c(3:8,10))
{
  df_normal[,n] <- as.numeric(as.character(df_normal[,n]))
}

openxlsx::write.xlsx(df_normal, file = "20220401_SummaryStatistics.xlsx", overwrite = T)


# Contraste complejo ------------------------------------------------------

constraste3 <- c("BENEFIT3_38","BOR_41", "PatronPD_123_47b")
table(x$BENEFIT3_38)
table(x$BOR_41)
table(x$PatronPD_123_47b)

normal_df2 <- NULL
dfvbase <- data.frame(matrix("-", ncol = 4, byrow = T))
names(dfvbase) <- c("Contraste","mean","sd","var")
dfv_stats <- NULL
dfvtuk <- data.frame(matrix("-", ncol = 7, byrow = T))
names(dfvtuk) <- c("diff","lwr","upr","p.adj", "group", "Variable","Contraste")

# http://www.sthda.com/english/wiki/one-way-anova-test-in-r

for(h in new_fact)
{
  # h <- new_fact[1]
  
  for(co in constraste3)
  {
    # co <- constraste[4]
    print(paste0("Contrast ", co, " - Variable ", h))
    temp_df <- x[,c(co, h)]
    temp2 <- temp_df[complete.cases(temp_df),]
    
    if(nrow(temp2) > 5 & table(temp2[,co])[1] > 1 & table(temp2[,co])[2] > 1 & var(temp2[,h]) != 0)
    {
      print(paste0("Contrast ", co, " - Variable ", h, "-----------------> SI"))
      
      table(temp2[,co])
      
      # library(dplyr)
      # temp2 %>% group_by_at (co) %>%
      #   summarise(
      #     count = n(),
      #     mean = mean(h, na.rm = TRUE),
      #     sd = sd(h, na.rm = TRUE)
      #   )
      
      dfv <- data.frame(aggregate(temp2[,h], list(temp2[,co]), mean),aggregate(temp2[,h], by =list(temp2[,co]), sd))
      names(dfv) <- c("Contraste", "mean", "ken","sd")
      dfv$var <- h
      dfv <- dfv[,-3]
      dfv$Contraste <- paste0(co, "_", dfv$Contraste)
      dfvbase <- rbind(dfvbase, dfv)
      
      temp3 <- temp2
      names(temp3) <- c("Contraste","Variable")
      
      library("ggpubr")
      ggboxplot(temp3, x = "Contraste", y = "Variable", 
                color = "Contraste", palette = c("#00AFBB", "#E7B800", "#FC4E07"),
                # order = c("ctrl", "trt1", "trt2"),
                ylab = h, xlab = "Contrast")
      
      library("ggpubr")
      ggline(temp3, x = "Contraste", y = "Variable", 
             add = c("mean_se", "jitter"), 
             # order = c("ctrl", "trt1", "trt2"),
             ylab = "Variable", xlab = "Contrast")
      
      
      
      # # Box plot
      # boxplot(Variable ~ Contraste, data = temp3,
      #         xlab = "Treatment", ylab = "Weight",
      #         frame = FALSE, col = c("#00AFBB", "#E7B800", "#FC4E07"))
      # # plotmeans
      # library("gplots")
      # plotmeans(Variable ~ Contraste, data = temp3, frame = FALSE,
      #           xlab = "Treatment", ylab = "Weight",
      #           main="Mean Plot with 95% CI") 
      
      
      # Compute the analysis of variance
      temp3$Contraste <- as.factor(temp3$Contraste)
      res.aov <- aov(Variable ~ Contraste, data = temp3)
      # Summary of the analysis
      s1 <- summary(res.aov)
      dfv_stats <- c(dfv_stats, h, co, unlist(s1[[1]][1,][5]))
      
      # Si anova es signifivcativo podemos mirar las diferencias entre grupos
      tuk <- TukeyHSD(res.aov)
      tuko <- data.frame(tuk$Contraste)
      tuko$group <- rownames(tuko)
      tuko$Variable <- h
      tuko$Contraste <- co
      
      
      dfvtuk <- rbind(dfvtuk, tuko)
      
      # p balio ajustatua
      pairwise.t.test(temp3$Variable, temp3$Contraste,
                      p.adjust.method = "BH")
      
      # 1. Homogeneity of variances
      plot(res.aov, 1)
      
      # We recommend Levene’s test, which is less sensitive to departures from normal distribution. 
      # The function leveneTest() [in car package] will be used:
      
      library(car)
      leveneTest(Variable ~ Contraste, data = temp3)
      # menos de 0.05 no hay homogeneidad
      
      
    }
    
  }
}

df_complejo <- data.frame(matrix(dfv_stats, ncol = 3, byrow = TRUE))
names(df_complejo) <- c("Variable","Contrast","Anova")

dfvbase <- dfvbase[-1,]
df_complejo <- df_complejo[-1,]
dfvtuk <- dfvtuk[-1, ]

dfvbase$mean <- as.numeric(as.character(dfvbase$mean))
dfvbase$sd <- as.numeric(as.character(dfvbase$sd))
df_complejo$Anova <- as.numeric(as.character(df_complejo$Anova))
dfvtuk$diff <- as.numeric(as.character(dfvtuk$diff))
dfvtuk$lwr <- as.numeric(as.character(dfvtuk$lwr))
dfvtuk$upr <- as.numeric(as.character(dfvtuk$upr))
dfvtuk$'p.adj' <- as.numeric(as.character(dfvtuk$'p.adj'))

dfvtuk <- dfvtuk[,c(7,5,6,1:4)]

openxlsx::write.xlsx(dfvbase, file = "20220401_Valores_var_complejas.xlsx", overwrite = T)
openxlsx::write.xlsx(df_complejo, file = "20220401_anova_var_complejas.xlsx", overwrite = T)
openxlsx::write.xlsx(dfvtuk, file = "20220401_tuk_var_complejas.xlsx", overwrite = T)




# Survival ----------------------------------------------------------------

x$FechaSIRT_2 <- as.Date(x$FechaSIRT_40, origin = "1899-12-30")
x$FechaPD_2 <- as.Date(x$FechaPD_46, origin = "1899-12-30")
x$FechaPFS_2 <- as.Date(x$FechaPFS_50, origin = "1899-12-30")
x$FechaOS_2 <- as.Date(x$FechaOS_52, origin = "1899-12-30")

grep("Fecha", names(x), value = T)
grep("OS_", names(x), value = T)
grep("PD_", names(x), value = T)
grep("PFS", names(x), value = T)

x$OS_days <- x$FechaOS_2 - x$FechaSIRT_2
x$PD_days <- x$FechaPD_2 - x$FechaSIRT_2
x$PDF_days <- x$FechaPFS_2 - x$FechaSIRT_2

x$OS_days <- x$OS_days_new_691
x$PD_days <- x$PD_days_new_692
x$PDF_days <- x$PFS_days_new_693

x$Muerte_51 <- as.numeric(as.character(x$Muerte_51))

library(foreign)
library(dplyr)
library(survival)
library(survminer)
library(ggpubr)

grep("Plasma",names(x))


# OS ----------------------------------------------------------------------

# Función para sacar los valores que tenemos de 0 y 1
taula_zenbat <- function(taula)
{
  izen <- names(taula)
  if(length(izen) == 2)
  {
    erantzuna <- unname(taula)
  } else if(length(izen) == 1)
  {
    if(izen == "0")
    {
      erantzuna <- c(unname(taula), 0)
    } else if(izen == "1")
    {
      erantzuna <- c(0, unname(taula))
    }
  }
  return(erantzuna)
}

survival_df <- NULL
j <- 0
for(h in new_fact)
{
  j <- j + 1
  print(paste0(j, ": ", h))
  
  # if(sum(is.na(x[,h])) < 40 & table(x[,h])[1] > 2 & !is.na((table(x[,h])[2] > 2 | !is.na(table(x[,h])[2]))))
  if(sum(is.na(x[,h])) < 35 & var(x[,h], na.rm = T))
  {
    x$hmean2 <-  x$hmedian2 <- x$hq1_234 <- x$hq123_4 <- hmean <- hmedian <- q <- NULL
    hmean <- mean(x[,h], na.rm = T)
    hmedian <- median(x[,h], na.rm = T)
    q <- quantile (x[,h], na.rm = T)
    x$hmean2 <- as.factor(ifelse(x[,h] <= hmean, 0, 1))
    x$hmedian2 <- as.factor(ifelse(x[,h] <= hmedian, 0, 1))
    x$hq1_234 <- as.factor(ifelse(x[,h] <= q[2], 0, 1))
    x$hq123_4 <- as.factor(ifelse(x[,h] <= q[4], 0, 1))
    
    
    sel_var <- c("OS_days","Muerte_51","hmean2")
    form <- as.formula("Surv(OS_days, Muerte_51) ~ hmean2")
    su <- coxph(form, data = x[complete.cases(x[,sel_var]),sel_var])
    ss <- summary(su)
    smean <- ss$coefficients[5]
    
    ggsurvplot(x,
               fit = survfit(Surv(OS_days, Muerte_51) ~ hmean2, data = x),
               xlab = "Days",
               ylab = "Overall survival probability",
               conf.int = FALSE,
               pval = TRUE)
    
    sel_var <- c("OS_days","Muerte_51","hmedian2")
    form <- as.formula("Surv(OS_days, Muerte_51) ~ hmedian2")
    su <- coxph(form, data = x[complete.cases(x[,sel_var]),sel_var])
    ss <- summary(su)
    smedian <- ss$coefficients[5]
    
    sel_var <- c("OS_days","Muerte_51","hq1_234")
    form <- as.formula("Surv(OS_days, Muerte_51) ~ hq1_234")
    su <- coxph(form, data = x[complete.cases(x[,sel_var]),sel_var])
    ss <- summary(su)
    s1_234 <- ss$coefficients[5]
    
    sel_var <- c("OS_days","Muerte_51","hq123_4")
    form <- as.formula("Surv(OS_days, Muerte_51) ~ hq123_4")
    su <- coxph(form, data = x[complete.cases(x[,sel_var]),sel_var])
    ss <- summary(su)
    s123_4 <- ss$coefficients[5]
    
    tau1 <- taula_zenbat(table(x$hmean2))
    tau2 <- taula_zenbat(table(x$hmedian2))
    tau3 <- taula_zenbat(table(x$hq1_234))
    tau4 <- taula_zenbat(table(x$hq123_4))
    
    survival_df <- c(survival_df, h, smean, tau1, smedian, tau2, s1_234, tau3, s123_4, tau4)
  }
}

df_survival <- data.frame(matrix(survival_df, ncol = 13, byrow = TRUE))
names(df_survival) <- c("Variable","Mean","Mean0","Mean1","Median","Med0","Med1",
                        "Q1_234","QA0","QA1","Q123_4","QB0","QB1")
head(df_survival)
for(n in 2:5)
{
  df_survival[,n] <- as.numeric(as.character(df_survival[,n]))
}

openxlsx::write.xlsx(df_survival, file = "20220401_SurvivalStatistics_OS.xlsx", overwrite = T)

# PD ----------------------------------------------------------------------

survival_df <- NULL
j <- 0
for(h in new_fact)
{
  j <- j + 1
  print(paste0(j, ": ", h))
  
  if(sum(is.na(x[,h])) < 35 & var(x[,h], na.rm = T))
  {
    x$hmean2 <-  x$hmedian2 <- x$hq1_234 <- x$hq123_4 <- hmean <- hmedian <- q <- NULL
    hmean <- mean(x[,h], na.rm = T)
    hmedian <- median(x[,h], na.rm = T)
    q <- quantile (x[,h], na.rm = T)
    x$hmean2 <- as.factor(ifelse(x[,h] <= hmean, 0, 1))
    x$hmedian2 <- as.factor(ifelse(x[,h] <= hmedian, 0, 1))
    x$hq1_234 <- as.factor(ifelse(x[,h] <= q[2], 0, 1))
    x$hq123_4 <- as.factor(ifelse(x[,h] <= q[4], 0, 1))
    
    sel_var <- c("PD_days","Muerte_51","hmean2")
    form <- as.formula("Surv(PD_days, Muerte_51) ~ hmean2")
    su <- coxph(form, data = x[complete.cases(x[,sel_var]),sel_var])
    ss <- summary(su)
    smean <- ss$coefficients[5]
    
    # ggsurvplot(x,
    #            fit = survfit(Surv(OS_days, Muerte_51) ~ hmean2, data = x),
    #            xlab = "Days",
    #            ylab = "Overall survival probability",
    #            conf.int = FALSE,
    #            pval = TRUE)
    
    sel_var <- c("PD_days","Muerte_51","hmedian2")
    form <- as.formula("Surv(PD_days, Muerte_51) ~ hmedian2")
    su <- coxph(form, data = x[complete.cases(x[,sel_var]),sel_var])
    ss <- summary(su)
    smedian <- ss$coefficients[5]
    
    sel_var <- c("PD_days","Muerte_51","hq1_234")
    form <- as.formula("Surv(PD_days, Muerte_51) ~ hq1_234")
    su <- coxph(form, data = x[complete.cases(x[,sel_var]),sel_var])
    ss <- summary(su)
    s1_234 <- ss$coefficients[5]
    
    sel_var <- c("PD_days","Muerte_51","hq123_4")
    form <- as.formula("Surv(PD_days, Muerte_51) ~ hq123_4")
    su <- coxph(form, data = x[complete.cases(x[,sel_var]),sel_var])
    ss <- summary(su)
    s123_4 <- ss$coefficients[5]
    
    survival_df <- c(survival_df, h, smean, smedian, s1_234, s123_4)
  }
}

df_survival <- data.frame(matrix(survival_df, ncol = 5, byrow = TRUE))
names(df_survival) <- c("Variable","Mean","Median","Q1_234","Q123_4")

head(df_survival)
for(n in 2:5)
{
  df_survival[,n] <- as.numeric(as.character(df_survival[,n]))
}

openxlsx::write.xlsx(df_survival, file = "20220401_SurvivalStatistics_PD.xlsx", overwrite = T)

# PDF ----------------------------------------------------------------------

survival_df <- NULL
j <- 0
for(h in new_fact)
{
  j <- j + 1
  print(paste0(j, ": ", h))
  
  if(sum(is.na(x[,h])) < 35 & var(x[,h], na.rm = T))
  {
    x$hmean2 <-  x$hmedian2 <- x$hq1_234 <- x$hq123_4 <- hmean <- hmedian <- q <- NULL
    hmean <- mean(x[,h], na.rm = T)
    hmedian <- median(x[,h], na.rm = T)
    q <- quantile (x[,h], na.rm = T)
    x$hmean2 <- as.factor(ifelse(x[,h] <= hmean, 0, 1))
    x$hmedian2 <- as.factor(ifelse(x[,h] <= hmedian, 0, 1))
    x$hq1_234 <- as.factor(ifelse(x[,h] <= q[2], 0, 1))
    x$hq123_4 <- as.factor(ifelse(x[,h] <= q[4], 0, 1))
    
    sel_var <- c("PDF_days","Muerte_51","hmean2")
    form <- as.formula("Surv(PDF_days, Muerte_51) ~ hmean2")
    su <- coxph(form, data = x[complete.cases(x[,sel_var]),sel_var])
    ss <- summary(su)
    smean <- ss$coefficients[5]
    
    
    
    
    # ggsurvplot(x,
    #            fit = survfit(Surv(OS_days, Muerte_51) ~ hmean2, data = x),
    #            xlab = "Days",
    #            ylab = "Overall survival probability",
    #            conf.int = FALSE,
    #            pval = TRUE)
    
    sel_var <- c("PDF_days","Muerte_51","hmedian2")
    form <- as.formula("Surv(PDF_days, Muerte_51) ~ hmedian2")
    su <- coxph(form, data = x[complete.cases(x[,sel_var]),sel_var])
    ss <- summary(su)
    smedian <- ss$coefficients[5]
    
    sel_var <- c("PDF_days","Muerte_51","hq1_234")
    form <- as.formula("Surv(PDF_days, Muerte_51) ~ hq1_234")
    su <- coxph(form, data = x[complete.cases(x[,sel_var]),sel_var])
    ss <- summary(su)
    s1_234 <- ss$coefficients[5]
    
    sel_var <- c("PDF_days","Muerte_51","hq123_4")
    form <- as.formula("Surv(PDF_days, Muerte_51) ~ hq123_4")
    su <- coxph(form, data = x[complete.cases(x[,sel_var]),sel_var])
    ss <- summary(su)
    s123_4 <- ss$coefficients[5]
    
    survival_df <- c(survival_df, h, smean, smedian, s1_234, s123_4)
  }
}

df_survival <- data.frame(matrix(survival_df, ncol = 5, byrow = TRUE))
names(df_survival) <- c("Variable","Mean","Median","Q1_234","Q123_4")

head(df_survival)
for(n in 2:5)
{
  df_survival[,n] <- as.numeric(as.character(df_survival[,n]))
}

openxlsx::write.xlsx(df_survival, file = "20220401_SurvivalStatistics_PDF.xlsx", overwrite = T)


# Plots -------------------------------------------------------------------



library(survminer)

hist(dat$Plasma.cells, breaks = 25)
dat$sel <- ifelse(dat$Plasma.cells < 0.001, 0, 1)
table(dat$sel)

ggsurvplot(dat,
           fit = survfit(Surv(OS_days, Muerte_51) ~ sel, data = dat),
           xlab = "Days",
           ylab = "Overall survival probability")

hist(dat$Only50RNA2, breaks = 25)
dat$sel <- ifelse(dat$Only50RNA2 < 20, 0, 1)
table(dat$sel)

form <- as.formula("Surv(OS_days, Muerte_51) ~ hmedian2")
su <- coxph(form, data = x)
ss <- summary(su)
smedian <- ss$coefficients[5]

ggsurvplot(x,
           fit = survfit(Surv(OS_days, Muerte_51) ~ hmedian2, data = x),
           xlab = "Days",
           ylab = "Overall survival probability",
           conf.int = TRUE,
           pval = TRUE)


# OS with BENEFIT2 --------------------------------------------------------


# Función para sacar los valores que tenemos de 0 y 1
taula_zenbat <- function(taula)
{
  izen <- names(taula)
  if(length(izen) == 2)
  {
    erantzuna <- unname(taula)
  } else if(length(izen) == 1)
  {
    if(izen == "0")
    {
      erantzuna <- c(unname(taula), 0)
    } else if(izen == "1")
    {
      erantzuna <- c(0, unname(taula))
    }
  }
  return(erantzuna)
}

smean_fin <- data.frame(matrix(NA, ncol = 5, byrow = T))
names(smean_fin) <- c("Coef","P-value","Group","Variable","Metric")
smedian_fin <- data.frame(matrix(NA, ncol = 5, byrow = T))
names(smedian_fin) <- c("Coef","P-value","Group","Variable","Metric")
s1_234_fin <- data.frame(matrix(NA, ncol = 5, byrow = T))
names(s1_234_fin) <- c("Coef","P-value","Group","Variable","Metric")
s123_4_fin <- data.frame(matrix(NA, ncol = 5, byrow = T))
names(s123_4_fin) <- c("Coef","P-value","Group","Variable","Metric")



survival_df <- NULL
j <- 0
for(h in new_fact)
{
  j <- j + 1
  print(paste0(j, ": ", h))
  
  # if(sum(is.na(x[,h])) < 40 & table(x[,h])[1] > 2 & !is.na((table(x[,h])[2] > 2 | !is.na(table(x[,h])[2]))))
  if(sum(is.na(x[,h])) < 35 & var(x[,h], na.rm = T))
  {
    x$hmean2 <-  x$hmedian2 <- x$hq1_234 <- x$hq123_4 <- hmean <- hmedian <- q <- NULL
    hmean <- mean(x[,h], na.rm = T)
    hmedian <- median(x[,h], na.rm = T)
    q <- quantile (x[,h], na.rm = T)
    x$hmean2 <- as.factor(ifelse(x[,h] <= hmean, 0, 1))
    x$hmedian2 <- as.factor(ifelse(x[,h] <= hmedian, 0, 1))
    x$hq1_234 <- as.factor(ifelse(x[,h] <= q[2], 0, 1))
    x$hq123_4 <- as.factor(ifelse(x[,h] <= q[4], 0, 1))
    
    
    sel_var <- c("OS_days","Muerte_51","hmean2","BENEFIT2_39")
    form <- as.formula("Surv(OS_days, Muerte_51) ~ hmean2 + BENEFIT2_39")
    su <- coxph(form, data = x[complete.cases(x[,sel_var]),sel_var])
    ss <- summary(su)
    smean <- ss$coefficients[5]
    smean <- data.frame(ss$coefficients[c(1,2),c(1,5)])
    smean$Group <- paste0(h, "_", rownames(smean)[2])
    smean$Variable <- c(h, rownames(smean)[2])
    names(smean)[1:2] <- c("Coef","P-value")
    smean$Metric <- "Mean"
    smean_fin <- rbind(smean_fin, smean)
    
    
    sel_var <- c("OS_days","Muerte_51","hmedian2","BENEFIT2_39")
    form <- as.formula("Surv(OS_days, Muerte_51) ~ hmedian2 + BENEFIT2_39")
    su <- coxph(form, data = x[complete.cases(x[,sel_var]),sel_var])
    ss <- summary(su)
    smedian <- data.frame(ss$coefficients[c(1,2),c(1,5)])
    smedian$Group <- paste0(h, "_", rownames(smedian)[2])
    smedian$Variable <- c(h, rownames(smedian)[2])
    names(smedian)[1:2] <- c("Coef","P-value")
    smedian$Metric <- "Median"
    smedian_fin <- rbind(smedian_fin, smedian)
    
    sel_var <- c("OS_days","Muerte_51","hq1_234","BENEFIT2_39")
    form <- as.formula("Surv(OS_days, Muerte_51) ~ hq1_234 + BENEFIT2_39")
    su <- coxph(form, data = x[complete.cases(x[,sel_var]),sel_var])
    ss <- summary(su)
    # s1_234 <- ss$coefficients[5]
    s1_234 <- data.frame(ss$coefficients[c(1,2),c(1,5)])
    s1_234$Group <- paste0(h, "_", rownames(s1_234)[2])
    s1_234$Variable <- c(h, rownames(s1_234)[2])
    names(s1_234)[1:2] <- c("Coef","P-value")
    s1_234$Metric <- "s1_234"
    s1_234_fin <- rbind(s1_234_fin, s1_234)
    
    sel_var <- c("OS_days","Muerte_51","hq123_4","BENEFIT2_39")
    form <- as.formula("Surv(OS_days, Muerte_51) ~ hq123_4 + BENEFIT2_39")
    su <- coxph(form, data = x[complete.cases(x[,sel_var]),sel_var])
    ss <- summary(su)
    # s123_4 <- ss$coefficients[5]
    s123_4 <- data.frame(ss$coefficients[c(1,2),c(1,5)])
    s123_4$Group <- paste0(h, "_", rownames(s123_4)[2])
    s123_4$Variable <- c(h, rownames(s123_4)[2])
    names(s123_4)[1:2] <- c("Coef","P-value")
    s123_4$Metric <- "s123_4"
    s123_4_fin <- rbind(s123_4_fin, s123_4)
    
    # tau1 <- taula_zenbat(table(x$hmean2))
    # tau2 <- taula_zenbat(table(x$hmedian2))
    # tau3 <- taula_zenbat(table(x$hq1_234))
    # tau4 <- taula_zenbat(table(x$hq123_4))
    # 
    # survival_df <- c(survival_df, h, smean, tau1, smedian, tau2, s1_234, tau3, s123_4, tau4)
  }
}

# df_survival <- data.frame(matrix(survival_df, ncol = 13, byrow = TRUE))
# names(df_survival) <- c("Variable","Mean","Mean0","Mean1","Median","Med0","Med1",
#                         "Q1_234","QA0","QA1","Q123_4","QB0","QB1")

ben_var <- rbind(smean_fin,smedian_fin, s1_234_fin, s123_4_fin)
ben_var <- ben_var[-1,]
ben_var <- ben_var[,c(3,4,5,1,2)]
ben_var$Coef <- as.numeric(as.character(ben_var$Coef))
ben_var$`P-value` <- as.numeric(as.character(ben_var$`P-value`))


openxlsx::write.xlsx(ben_var, file = "20220401_SurvivalStatistics_OS_with_BENEFIT2.xlsx", overwrite = T)


# cox with more covariates ------------------------------------------------

smean_fin <- data.frame(matrix(NA, ncol = 5, byrow = T))
names(smean_fin) <- c("Coef","P-value","Group","Variable","Metric")
smedian_fin <- data.frame(matrix(NA, ncol = 5, byrow = T))
names(smedian_fin) <- c("Coef","P-value","Group","Variable","Metric")
s1_234_fin <- data.frame(matrix(NA, ncol = 5, byrow = T))
names(s1_234_fin) <- c("Coef","P-value","Group","Variable","Metric")
s123_4_fin <- data.frame(matrix(NA, ncol = 5, byrow = T))
names(s123_4_fin) <- c("Coef","P-value","Group","Variable","Metric")

multicontraste <- c("BENEFIT2_39","PVI_29","TACE_30", "Etiology_32", "AgeMay65_33" , "AFPMzy400_35")

extract_surv <- function(sel_var, x, concept)
{
  # sel_var = dos variables de supervivencia y el resto
  # x = dataset con los datos
  
  form <- as.formula(paste0("Surv(", sel_var[1], ", ", sel_var[2], ") ~ ", paste0(sel_var[-c(1,2)], collapse = " + ")))
  su <- coxph(form, data = x[complete.cases(x[,sel_var]),sel_var])
  ss <- summary(su)
  # smean <- ss$coefficients[5]
  scolen <- nrow(ss$coefficients)
  smean <- data.frame(ss$coefficients[c(1:scolen),c(1,5)])
  smean$Group <- paste(sel_var[-c(1,2)], collapse = "_")
  smean$Variable <- rownames(smean)
  names(smean)[1:2] <- c("Coef","P-value")
  smean$Metric <- concept
  return(smean)
}

# all posibilities
posib <- expand.grid(0:5, 0:5, 0:5, 0:5, 0:5)
posib <- posib[-1,]
posib$summary <- 0

for(aa in 1:nrow(posib))
{
  bb <- sort(unique(unlist(posib[aa, ])))
  cc <- paste(bb,collapse="")
  posib$summary[aa] <- gsub("0","",cc)
}
un_pos <- unique(posib$summary)

# para saber qué covariates entran en el análisis
cual_cov <- function(vector, string)
{
  s <- 1
  vec <- NULL
  for(h in 1:nchar(string))
  {
    vec_temp <- as.numeric(substr(string,s, s))
    vec <- c(vec, vec_temp)
    s <- s + 1
  }
  vecnname <- vector[vec]
  return(vecnname)
}


survival_df <- NULL
j <- 0
for(h in new_fact[-c(1:4)])
{
  j <- j + 1
  print(paste0(j, ": ", h))
  
  # if(sum(is.na(x[,h])) < 40 & table(x[,h])[1] > 2 & !is.na((table(x[,h])[2] > 2 | !is.na(table(x[,h])[2]))))
  if(sum(is.na(x[,h])) < 25 & var(x[,h], na.rm = T))
  {
    x$hmean2 <-  x$hmedian2 <- x$hq1_234 <- x$hq123_4 <- hmean <- hmedian <- q <- NULL
    hmean <- mean(x[,h], na.rm = T)
    hmedian <- median(x[,h], na.rm = T)
    q <- quantile (x[,h], na.rm = T)
    x$hmean2 <- as.factor(ifelse(x[,h] <= hmean, 0, 1))
    x$hmedian2 <- as.factor(ifelse(x[,h] <= hmedian, 0, 1))
    x$hq1_234 <- as.factor(ifelse(x[,h] <= q[2], 0, 1))
    x$hq123_4 <- as.factor(ifelse(x[,h] <= q[4], 0, 1))
    
    for(coco in un_pos)
    {
      coco_var <- cual_cov(vector = multicontraste, string = coco)
      sel_var <- c("OS_days","Muerte_51","hmean2",coco_var)
      smean <- extract_surv(sel_var, x, concept = paste0("Mean_", h))
      smean_fin <- rbind(smean_fin, smean)
      
      sel_var <- c("OS_days","Muerte_51","hmedian2",coco_var)
      smedian <- extract_surv(sel_var, x, concept = paste0("Median", h))
      smedian_fin <- rbind(smedian_fin, smedian)
      
      sel_var <- c("OS_days","Muerte_51","hq1_234",coco_var)
      s1_234 <- extract_surv(sel_var, x, concept = paste0("s1_234", h))
      s1_234_fin <- rbind(s1_234_fin, s1_234)
      
      sel_var <- c("OS_days","Muerte_51","hq123_4",coco_var)
      s123_4 <- extract_surv(sel_var, x, concept = paste0("s123_4", h))
      s123_4_fin <- rbind(s123_4_fin, s123_4)
    }
    
    
    
  }
}

# df_survival <- data.frame(matrix(survival_df, ncol = 13, byrow = TRUE))
# names(df_survival) <- c("Variable","Mean","Mean0","Mean1","Median","Med0","Med1",
#                         "Q1_234","QA0","QA1","Q123_4","QB0","QB1")

ben_var <- rbind(smean_fin,smedian_fin, s1_234_fin, s123_4_fin)
ben_var <- ben_var[-1,]
ben_var <- ben_var[,c(3,4,5,1,2)]
ben_var$Coef <- as.numeric(as.character(ben_var$Coef))
ben_var$`P-value` <- as.numeric(as.character(ben_var$`P-value`))


openxlsx::write.xlsx(ben_var, file = "20220401_SurvivalStatistics_OS_with_combinations.xlsx", overwrite = T)


# EJemplo individual ------------------------------------------------------

h <- "P3_CK_IL_2_102"
library(survival)

j <- j + 1
print(paste0(j, ": ", h))

x$hmean2 <-  x$hmedian2 <- x$hq1_234 <- x$hq123_4 <- hmean <- hmedian <- q <- NULL
hmean <- mean(x[,h], na.rm = T)
hmedian <- median(x[,h], na.rm = T)
q <- quantile (x[,h], na.rm = T)
x$hmean2 <- as.factor(ifelse(x[,h] <= hmean, 0, 1))
x$hmedian2 <- as.factor(ifelse(x[,h] <= hmedian, 0, 1))
x$hq1_234 <- as.factor(ifelse(x[,h] <= q[2], 0, 1))
x$hq123_4 <- as.factor(ifelse(x[,h] <= q[4], 0, 1))

form <- as.formula("Surv(OS_days, Muerte_51) ~ hmedian2")
sel_var <- c("OS_days","Muerte_51","hmedian2","BENEFIT2_39")
table(x$hmedian2, x$BENEFIT2_39)
datu <- x[complete.cases(x[,sel_var]),sel_var]
su <- coxph(form, data = datu)
ss <- summary(su)
# smean <- ss$coefficients[5]

########################

sel_var <- c("OS_days","Muerte_51","hmean2", "BENEFIT2_39")
form <- as.formula("Surv(OS_days, Muerte_51) ~ hmean2 + BENEFIT2_39")
su <- coxph(form, data = x[complete.cases(x[,sel_var]),sel_var])
ss <- summary(su)
smean <- ss$coefficients[5]

ggsurvplot(x,
           fit = survfit(Surv(OS_days, Muerte_51) ~ hmean2 + BENEFIT2_39, data = x),
           xlab = "Days",
           ylab = "Overall survival probability",
           conf.int = FALSE,
           pval = TRUE)

############################################################

sel_var <- c("OS_days","Muerte_51","hmedian2", "BENEFIT2_39")
form <- as.formula("Surv(OS_days, Muerte_51) ~ hmedian2")
su <- coxph(form, data = x[complete.cases(x[,sel_var]),sel_var])
ss <- summary(su)
smean <- ss$coefficients[5]

ggsurvplot(x,
           fit = survfit(Surv(OS_days, Muerte_51) ~ hmedian2, data = x),
           xlab = "Days",
           ylab = "Overall survival probability",
           conf.int = FALSE,
           pval = TRUE)

############################################################

sel_var <- c("OS_days","Muerte_51","hmedian2", "BENEFIT2_39")
form <- as.formula("Surv(OS_days, Muerte_51) ~ BENEFIT2_39")
su <- coxph(form, data = x[complete.cases(x[,sel_var]),sel_var])
ss <- summary(su)
smean <- ss$coefficients[5]

ggsurvplot(x,
           fit = survfit(Surv(OS_days, Muerte_51) ~ BENEFIT2_39, data = x),
           xlab = "Days",
           ylab = "Overall survival probability",
           conf.int = FALSE,
           pval = TRUE)

############################################################

sel_var <- c("OS_days","Muerte_51","hmedian2", "BENEFIT2_39", "PVI_29", "TACE_30","AgeMay65_33")
form <- as.formula("Surv(OS_days, Muerte_51) ~ hq1_234 + BENEFIT2_39 + PVI_29 + TACE_30 + AgeMay65_33")
su <- coxph(form, data = x[complete.cases(x[,sel_var]),sel_var])
ss <- summary(su)
smean <- ss$coefficients[5]

ggsurvplot(x,
           fit = survfit(Surv(OS_days, Muerte_51) ~ hq1_234 + BENEFIT2_39 + PVI_29 + TACE_30 + AgeMay65_33, data = x),
           xlab = "Days",
           ylab = "Overall survival probability",
           conf.int = FALSE,
           pval = TRUE)

