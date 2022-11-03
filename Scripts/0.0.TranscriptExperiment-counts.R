
# Phaeobacter 


# Import data
CFU.abundance <- readxl::read_excel("Data/newcounttranscripts20220622.xlsx", sheet = 3)
CFU.abundance <- readxl::read_excel("Data/TranscriptCFU.xlsx", sheet = 1)

#CFU data
sum = CFU.abundance %>%
  group_by(day, count_of, culture, condition) %>%
  dplyr::summarize(count = n(),
            mean.log = mean(log_cells),
            sd.log = sd(log_cells))


# Plots

# Phaeobacter growth 
system_color2 = c("#9C3206", "#00A091", "black")
sum$culture = plyr::revalue(sum$culture, c("c"="Control", "p"="TDA", "d" = "NoTDA"))
sum = dplyr::rename(sum, System = culture)
sum$System = factor(sum$System, levels= c("TDA", "NoTDA", "Control"))

ggplot(sum[sum$count_of != "c",], aes(day, mean.log, color = System)) +
  geom_line(size=0.5, alpha = 0.6, linetype = 'dashed') + 
  geom_point(size=2) + 
  geom_errorbar(aes(ymax = mean.log + sd.log, ymin = mean.log - sd.log),
                alpha = 0.6, width = 0, size = 0.5) + 
  labs(x="\n Time (days) \n", y= "log10(cells/mL) \n") + 
  theme_bw(base_size = 12) +
  facet_grid(System~condition)+
  theme(legend.text = element_text(size = 11), 
        legend.title = element_text(size=12, face="bold")) +
  # scale_color_manual(values=system_color2, labels = c(expression(italic("P. inhibens (WT)")),expression(Delta*"tdaB"), "Control"))+
  # scale_fill_manual(values=system_color2, labels = c(expression(italic("P. inhibens (WT)")),expression(Delta*"tdaB"), "Control"))+
  scale_x_continuous(breaks = c(0, 1,2,3,4, 7))+
  theme(axis.text.y = element_text(colour = "black", face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold"), 
        legend.text = element_text(colour ="black"),
        legend.position = "right", axis.title.y = element_text(face = "bold"), 
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title.x = element_text(face = "bold", colour = "black"), 
        legend.title = element_text(colour = "black", face = "bold"), 
        # panel.grid.major = element_blank(),
        # panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))

ggplot(sum[sum$count_of == "c",], aes(day, mean.log, color = System)) +
  geom_line(size=0.5, alpha = 0.6, linetype = 'dashed') + 
  geom_point(size=2) + 
  geom_errorbar(aes(ymax = mean.log + sd.log, ymin = mean.log - sd.log),
                alpha = 0.6, width = 0, size = 0.5) + 
  labs(x="\n Time (days) \n", y= "log10(cells/mL) \n") + 
  theme_bw(base_size = 12) +
  facet_grid(System~condition)+
  theme(legend.text = element_text(size = 11), 
        legend.title = element_text(size=12, face="bold")) +
  # scale_color_manual(values=system_color2, labels = c(expression(italic("P. inhibens (WT)")),expression(Delta*"tdaB"), "Control"))+
  # scale_fill_manual(values=system_color2, labels = c(expression(italic("P. inhibens (WT)")),expression(Delta*"tdaB"), "Control"))+
  scale_x_continuous(breaks = c(0, 1,2,3,4, 7))+
  theme(axis.text.y = element_text(colour = "black", face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold"), 
        legend.text = element_text(colour ="black"),
        legend.position = "right", axis.title.y = element_text(face = "bold"), 
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title.x = element_text(face = "bold", colour = "black"), 
        legend.title = element_text(colour = "black", face = "bold"), 
        # panel.grid.major = element_blank(),
        # panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))

#ggsave(filename = "Fig1C-Algae.pdf", plot = S1, width=15, height=10, units="cm", dpi = 400)

# Total community 
C.abundance <- readxl::read_excel("Data/newcounttranscripts20220622.xlsx", sheet = 5)
sum.c = C.abundance %>%
  group_by(culture, condition, day) %>%
  dplyr::summarize(count = n(),
                   mean.log = mean(log_TC),
                   sd.log = sd(log_TC))

ggplot(sum.c, aes(, mean.log, color = culture)) +
  #geom_line(size=0.5, position = position_dodge(width = 2.5), alpha = 0.6, linetype = 'dashed') + 
  geom_point(size=2, position = position_dodge(width = 2.5)) + 
  geom_errorbar(aes(ymax = mean.log + sd.log, ymin = mean.log - sd.log),
                alpha = 0.6, position = position_dodge(width = 2.5), width = 0, size = 0.5) + 
  labs(x="\n Time (days) \n", y= "log10(cells/mL) \n") + 
  theme_bw(base_size = 12) +
  facet_grid(.~condition)+
  theme(legend.text = element_text(size = 11), 
        legend.title = element_text(size=12, face="bold")) +
  # scale_color_manual(values=system_color2, labels = c(expression(italic("P. inhibens (WT)")),expression(Delta*"tdaB"), "Control"))+
  # scale_fill_manual(values=system_color2, labels = c(expression(italic("P. inhibens (WT)")),expression(Delta*"tdaB"), "Control"))+
  # scale_x_continuous(breaks = c(0, 7, 14, 21, 28, 35, 42, 49, 56, 63, 70))+
  theme(axis.text.y = element_text(colour = "black", face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold"), 
        legend.text = element_text(colour ="black"),
        legend.position = "right", axis.title.y = element_text(face = "bold"), 
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title.x = element_text(face = "bold", colour = "black"), 
        legend.title = element_text(colour = "black", face = "bold"), 
        # panel.grid.major = element_blank(),
        # panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))

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

# TS abundance

NH4.abundance <- readxl::read_excel("Data/newcounttranscripts20220622.xlsx", sheet = 1)
F2.abundance <- readxl::read_excel("Data/newcounttranscripts20220622.xlsx", sheet = 2)

#NH4
  
  str(NH4.abundance)
  NH4.abundance = dplyr::rename(NH4.abundance, System = culture)
  
  NH4.sum = NH4.abundance %>%
    group_by(System, day) %>%
    summarise(count = n(),
              mean.log = mean(log10),
              sd.log = sd(log10))
str(NH4.sum)  

# Plot 

TSAbun_transc = ggplot(NH4.sum, aes(day, mean.log, color = System, fill = System)) +
  geom_point(data = NH4.abundance, aes(x = day, y = log10, group = System),
             alpha = 0.6, position = position_dodge(width = 1), size = 1.5, stroke = 0)+
  geom_line(size=0.5, position = position_dodge(width = 1), alpha = 0.6, linetype = 'dashed') + 
  geom_point(size=3, position = position_dodge(width = 1)) + 
  geom_errorbar(aes(ymax = mean.log + sd.log, ymin = mean.log - sd.log),
                alpha = 0.6, position = position_dodge(width = 1), width = 0, size = 0.5) + 
  labs(x="\n Time (days) \n", y= "log10(cells/mL) \n") + 
  theme_bw(base_size = 12) +
  theme(legend.text = element_text(size = 11), 
        legend.title = element_text(size=12, face="bold")) +
  labs(x="\n Time (days) \n", y= "log10(cells/mL) \n") + 
  theme_bw(base_size = 12) +
  theme(legend.text = element_text(size = 11), 
        legend.title = element_text(size=12, face="bold")) +
  scale_color_manual(values=system_color2, labels = c(expression(italic("P. inhibens (WT)")),expression(Delta*"tdaB"), "Control"))+
  scale_fill_manual(values=system_color2, labels = c(expression(italic("P. inhibens (WT)")),expression(Delta*"tdaB"), "Control"))+
  scale_x_continuous(breaks = c(0, 1,4,7))+
  scale_y_continuous(breaks = c(4,4.5,5,5.5,6.0))+
  theme(axis.text.y = element_text(colour = "black", face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold"), 
        legend.text = element_text(colour ="black"),
        legend.position = c(0.75,0.2), axis.title.y = element_text(face = "bold"), 
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title.x = element_text(face = "bold", colour = "black"), 
        legend.title = element_text(colour = "black", face = "bold"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))

ggsave(filename = "TSAbun_transc.pdf", plot = TSAbun_transc, width=15, height=10, units="cm", dpi = 400)



  
  
  
#F2

str(F2.abundance)

F2.sum = F2.abundance %>%
  group_by(culture, day) %>%
  summarise(count = n(),
            mean.log = mean(log10),
            sd.log = sd(log10))
str(F2.sum)  

# Plot 

ggplot(F2.sum, aes(day, mean.log, color = culture, fill = culture)) +
  geom_line(size=0.5, alpha = 0.6) + 
  geom_point(size=3) + 
  geom_ribbon(aes(ymax = mean.log + sd.log, ymin = mean.log - sd.log),
              alpha = 0.2, color = NA) + 
  labs(x="\n Time (days) \n", y= "log10(cells/mL) \n") + 
  theme_bw(base_size = 12) +
  theme(legend.text = element_text(size = 11), 
        legend.title = element_text(size=12, face="bold")) +
  #scale_color_manual(values=system_color2, labels = c(expression(italic("P. inhibens (WT)")),expression(Delta*"tdaB"), "Control"))+
  #scale_fill_manual(values=system_color2, labels = c(expression(italic("P. inhibens (WT)")),expression(Delta*"tdaB"), "Control"))+
  scale_x_continuous(breaks = c(0,1,4,7))+
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
