
# Tetraselmis suecica abundance over time


# Import data
abundance <- read.csv("Data/Abundance_Tetraselmis.csv", header = T, sep = ";")

# Data wrangling

# Average techrep (sheet 1) #

algae_tech <- abundance %>%
  group_by(time, day, culture, biorep) %>%
  dplyr::summarize(techrep_count=n(),
                   cellsml_mean = mean(cellsml),
                   cellsml_std = sd(cellsml),
                   log_mean = mean(log10),
                   log_std = sd(log10))


# Start concentration 
algae_start = abundance %>%
  group_by(time) %>%
  dplyr::summarize(count = n(),
            mean.start.conc = mean(cellsml),
            sd.start.conc = sd(cellsml))

# Average biorep  (sheet 1)#

algae_bio <- algae_tech %>%
  group_by(time, day, culture) %>%
  dplyr::summarize(biorep_count = n(),
                   cellsml_bio_mean = mean(cellsml_mean),
                   cellsml_bio_std = sd(cellsml_mean),
                   log_bio_mean = mean(log_mean),
                   log_bio_std = sd(log_mean))

# Plots

# mean tech reps #

ggplot(algae_tech, aes(day, log_mean, color = biorep))+
  geom_point(size=1) +
  geom_line(size=0.7) + 
  geom_errorbar(aes(ymax = log_mean + log_std, ymin = log_mean - log_std), width = 1) + 
  facet_grid(culture~.) +
  labs(x="Time (days)", y= "log(cells/ml)") + 
  theme_bw() +
  scale_x_continuous(breaks = seq(0, 70, by=7))

#Mean bioreps cultures #
system_color2 = c("#9C3206", "#00A091", "black")
algae_bio$culture = plyr::revalue(algae_bio$culture, c("c"="Control", "p"="TDA", "d" = "NoTDA"))
algae_bio = dplyr::rename(algae_bio, System = culture)
algae_bio$System = factor(algae_bio$System, levels= c("TDA", "NoTDA", "Control"))

algae_tech$culture = plyr::revalue(algae_tech$culture, c("c"="Control", "p"="TDA", "d" = "NoTDA"))
algae_tech = dplyr::rename(algae_tech, System = culture)
algae_tech$System = factor(algae_tech$System, levels= c("TDA", "NoTDA", "Control"))


S1 = ggplot(algae_bio, aes(day, log_bio_mean, color = System)) +
  geom_point(data = algae_tech, aes(x = day, y = log_mean, group = System),
             alpha = 0.6, position = position_dodge(width = 2.5), size = 1.5, stroke = 0)+
  geom_line(size=0.5, position = position_dodge(width = 2.5), alpha = 0.6, linetype = 'dashed') + 
  geom_point(size=3, position = position_dodge(width = 2.5)) + 
  geom_errorbar(aes(ymax = log_bio_mean + log_bio_std, ymin = log_bio_mean - log_bio_std),
                alpha = 0.6, position = position_dodge(width = 2.5), width = 0, size = 0.5) + 
  labs(x="\n Time (days) \n", y= "log10(cells/mL) \n") + 
  theme_bw(base_size = 12) +
  theme(legend.text = element_text(size = 11), 
        legend.title = element_text(size=12, face="bold")) +
  scale_color_manual(values=system_color2, labels = c(expression(italic("P. inhibens (WT)")),expression(Delta*"tdaB"), "Control"))+
  scale_fill_manual(values=system_color2, labels = c(expression(italic("P. inhibens (WT)")),expression(Delta*"tdaB"), "Control"))+
  scale_x_continuous(breaks = c(0, 7, 14, 21, 28, 35, 42, 49, 56, 63, 70))+
  theme(axis.text.y = element_text(colour = "black", face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold"), 
        legend.text = element_text(colour ="black"),
        legend.position = c(0.8,0.2), axis.title.y = element_text(face = "bold"), 
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title.x = element_text(face = "bold", colour = "black"), 
        legend.title = element_text(colour = "black", face = "bold"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))
S1 
#ggsave(filename = "algae_abundance_culture.png", plot=start_bio_plot_N, device="png", height=5, width=10, units="in", dpi=500)

ggsave(filename = "Fig1C-Algae.pdf", plot = S1, width=15, height=10, units="cm", dpi = 400)

#############
### Stats ###
#############


# Linear mixed model (initial model with interaction)
lmm.ts0 <- lmer(cellsml_mean ~ as.factor(culture) *as.factor(day) +(1|biorep), data=algae_tech)

# Check studentized residuals for outliers
hist(rstudent(lmm.ts0))

# Find values with studentized residuals above 2.5 (above 2 is a reasonable cutoff)
highResIndx=which(abs(rstudent(lmm.ts0))>2.5)

# Linear mixed model with interaction and removed outliers
lmm.ts <- lmer(cellsml_mean ~ as.factor(culture) *as.factor(day) +(1|biorep), data=algae_tech[-highResIndx,])

# Should random effect be included
summary(lmm.ts) # lineage has an effect so should be a random effect (not zero)

#Test model
Anova(lmm.ts) #interaction is borderline, time is significant

#Post-hoc analysis
emm.ts = emmeans(lmm.ts, ~ as.factor(culture)| as.factor(day) )
contrast(emm.ts, interaction = "pairwise", adjust = "Bonferroni") #day 49 has significance vs TDA


 ##############################################################################################################################################

 ## Phaeobacter abundance #### 

 
#### Import data ####

#phaeo <- read_xlsx(path = "NC-LTEE1_CFUCount.xlsx", sheet = 2, col_names = TRUE)
#write.table(phaeo, file = "NC-LTEE1_CFUCount.csv", row.names = TRUE, sep = ";")

pCFU <- read.csv(file = "NC-LTEE1_CFUCount.csv", header = T, sep = ";")

# data without day15 (transfer numbers) #

pCFU2<- pCFU[ which(pCFU$day !='15'),]


#### Data wrangling ####

# Average techrep # 

pCFU_tech <- pCFU2 %>%
  group_by(lineage, day) %>%
  dplyr::summarize(techrep_count = n(),
                   logcfuml_mean = mean(logcfuml),
                   logcfuml_std = sd(logcfuml),
                   cfuml_mean = mean(cfuml),
                   cfuml_std = sd(cfuml))

# Average biorep #

pCFU_bio <- pCFU_tech %>%
  group_by(day) %>%
  dplyr::summarize(biorep_count = n(),
                   logcfuml_bio_mean = mean(logcfuml_mean),
                   logcfuml_bio_std = sd(logcfuml_mean),
                   cfuml_bio_mean = mean(cfuml_mean),
                   cfuml_bio_std = sd(cfuml_mean))


#### Plots log ####

# P culture (average of lineages) # 

pCFU_plot <- ggplot(pCFU_bio, aes(day, logcfuml_bio_mean)) +
  geom_point(size=1) + 
  geom_line(size=0.3) + 
  geom_errorbar(aes(ymax = logcfuml_bio_mean + logcfuml_bio_std, ymin = logcfuml_bio_mean - logcfuml_bio_std), width = 1) + 
  labs(title = "Abundance of P. inhibens", x="Time (days)", y= "log(CFU/ml)") + 
  scale_x_continuous(breaks = seq(0, 70, by=7))+
  theme_bw()+
  scale_color_manual(values=wes_palette(name="Cavalcanti1"))

pCFU_plot
ggsave(filename = "CFU_P.png", plot=pCFU_plot, device="png", height=5, width=8, units="in", dpi=500)



# Lineages #

pCFU_lin_plot <- ggplot(pCFU_tech,aes(day, logcfuml_mean, color = lineage)) +
  geom_point(size=1) +
  geom_line() +
  geom_errorbar(aes(ymax = logcfuml_mean + logcfuml_std, ymin = logcfuml_mean - logcfuml_std), width = 1) +
  labs(title = "Abundance of P. inhibens in lineages", subtitle = "Lineages", x="Time (days)", y= "log(CFU/ml)") + 
  facet_grid(lineage~.)+
  scale_x_continuous(breaks = seq(0, 70, by=7))+
  theme_bw()+
  scale_color_manual(values=wes_palette(n = 4, name="Cavalcanti1"))
  
pCFU_lin_plot

ggsave(filename = "CFU_P_lineages.png", plot=pCFU_lin_plot, device="png", height=5, width=8, units="in", dpi=500)


#### Plots ####

# P culture (average of lineages) # 

pCFU_plot_N <- ggplot(pCFU_bio, aes(day, cfuml_bio_mean)) +
  geom_point(size=1) + 
  geom_line(size=0.3) + 
  geom_errorbar(aes(ymax = cfuml_bio_mean + cfuml_bio_std, ymin = cfuml_bio_mean - cfuml_bio_std), width = 1) + 
  labs(title = "Abundance of P. inhibens", x="Time (days)", y= "CFU/ml") + 
  scale_x_continuous(breaks = seq(0, 70, by=7))+
  theme_bw()+
  scale_color_manual(values=wes_palette(n=1, name="GrandBudapest1"))

pCFU_plot_N


# Lineages #

pCFU_lin_plot_N <- ggplot(pCFU_tech,aes(day, cfuml_mean, color = lineage)) +
  geom_point(size=1) +
  geom_line(size=0.3) +
  geom_errorbar(aes(ymax = cfuml_mean + cfuml_std, ymin = cfuml_mean - cfuml_std), width = 1) +
  labs(title = "Abundance of P. inhibens in lineages", subtitle = "Lineages", x="Time (days)", y= "CFU/ml") + 
  facet_grid(lineage~.)+
  scale_x_continuous(breaks = seq(0, 70, by=7))+
  theme_bw()+
  scale_color_manual(values=wes_palette(n=4, name="GrandBudapest1"))

pCFU_lin_plot_N
  
  