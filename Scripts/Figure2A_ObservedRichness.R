##########################################################

                  ### Alphadiversity ###

#########################################################

#################
### Load data ###
#################

load("Data/ReadO/ps.asv.reduced.wTree.RData")


#########################
### Observed richness ###
#########################

#Make alpha estimation (Observed)
set.seed(41)
observed = estimate_richness(ps.new, measures = "Observed")
# Make df
# Get sample ID as rows
rownames(observed) = sample_names(ps.new) 

# Extract sample_data fra PS
s <- data.frame(sample_data(ps.new)) 

# Get a vector with SampleID - for merging
SampleID <- row.names(observed)

#combine the observed with the sampleID column
even_stats <- data.frame(SampleID, observed) 

# merge metadata and observed  by row.names (by = 0)
alphadiv <- merge(even_stats, s, by = 0) 

# Rename factor, order level, add day and lineage
alphadiv = alphadiv %>% dplyr::rename(System = exp.condition)
alphadiv$System = factor(alphadiv$System, levels= c("TDA", "NoTDA", "Control"))
alphadiv$day = factor(alphadiv$day, levels = c("0", "7", "14", "21", "28", "35", "42", "49", "56", "63", "70"))
alphadiv$lineage = paste(alphadiv$System,alphadiv$bio.rep, sep = "")

# Summarize
alphadiv_stats = alphadiv %>%
  group_by(day, System) %>%
  summarise(count = n(),
            observed.mean = mean(Observed),
            observed.std = sd(Observed))



table.observed = alphadiv_stats %>%
  arrange(System, day)

write.csv(table.observed, file = "Alphadiversity-table.csv")


##############
### Figure ###
##############

# color palette
system_color2 = c("#9C3206", "#00A091", "black")
# Figure
Fig3A = ggplot(alphadiv_stats , aes(x=day, y= observed.mean, color = System, fill = System))+
  geom_errorbar(aes(ymin=observed.mean-observed.std, ymax=observed.mean+observed.std, group = System), size = 0.5,width=0,
                position=position_dodge(0.5), alpha = 0.6)+
  geom_point(data = alphadiv, aes(x = day, y = Observed, group = System),
             alpha = 0.6, position = position_dodge(width = 0.5), size = 2, stroke = 0)+
  geom_point(aes(x = day, y = observed.mean, group = System),
             alpha = 1, position = position_dodge(width = 0.5), size = 3)+
  theme_bw(base_size = 10) +
  scale_color_manual(values=system_color2, labels = c(expression(italic("P. inhibens (WT)")),expression(Delta*"tdaB"), "Control"))+
  scale_fill_manual(values=system_color2, labels = c(expression(italic("P. inhibens (WT)")),expression(Delta*"tdaB"), "Control"))+
  labs(y = "Observed richness", x="Time (days)")+
  #theme(strip.text.x = element_blank()) +
  theme(legend.text = element_text(colour = "black"),
        legend.title = element_text(face = "bold",colour = "black"),
        legend.position= "right")+
  theme(axis.text = element_text(color = "Black", face = "bold"),
        axis.title = element_text(color = "Black", face = "bold"))+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
Fig3A 

legend.richness<- as_ggplot(ggpubr::get_legend(Fig3A))

ggsave("Fig3/Fig3A.pdf",
       plot = Fig3A,
       units="cm",
       width=17,
       height=7,
       dpi = 300)

ggsave("Fig3/Fig3A.lengend.pdf",
       plot = legend.richness,
       units="cm",
       width=4,
       height=3,
       dpi = 300)

#############
### Stats ###
#############

# Linear mixed model (initial model with interaction)
lmm.alpha0 <- lmer(Observed ~ System *day +(1|lineage), data=alphadiv)

# Check studentized residuals for outliers
hist(rstudent(lmm.alpha0))

# Find values with studentized residuals above 2.5 (above 2 is a reasonable cutoff)
highResIndx=which(abs(rstudent(lmm.alpha0))>2.5)

# Linear mixed model with interaction and removed outliers
lmm.alpha <- lmer(Observed ~ System *day +(1|lineage), data=alphadiv[-highResIndx,])

# Should random effect be included
summary(lmm.alpha) # lineage has an effect so should be a random effect (not zero)

#Test model
Anova(lmm.alpha) #interaction is borderline, time is significant

#Post-hoc analysis
emm.alpha = emmeans(lmm.alpha, ~ System| day )
contrast(emm.alpha, interaction = "pairwise", adjust = "Bonferroni") #day 42 has significance vs TDA

# Get measures 
sum.alphadiv = alphadiv %>%
  group_by(day, System) %>%
  summarise(mean.observed = mean(Observed),
            sd.observed = sd(Observed),
            count = n())

write.csv(sum.alphadiv, file = "observedrichness_table.csv")


# PD
PD <- phyloseq_phylo_div(ps.new, measures = c("PD"))
alpha.PD <- merge(PD, sample_data(ps.new), by = "row.names")

# Rename factor, order level, add day and lineage
alpha.PD = alpha.PD  %>% dplyr::rename(System = exp.condition)
alpha.PD$System = factor(alpha.PD$System, levels= c("TDA", "NoTDA", "Control"))
alpha.PD$day = factor(alpha.PD$day, levels = c("0", "7", "14", "21", "28", "35", "42", "49", "56", "63", "70"))
alpha.PD$lineage = paste(alpha.PD$System,alpha.PD$bio.rep, sep = "")

# Summarize
alpha.PD_stats = alpha.PD %>%
  group_by(day, System) %>%
  summarise(count = n(),
            PD.mean = mean(PD),
            PD.std = sd(PD))

table.PD = alpha.PD_stats %>%
  arrange(System, time)


##############
### Figure ###
##############

# color palette
system_color2 = c("#9C3206", "#00A091", "black")

# Figure
Fig3B = ggplot(alpha.PD_stats, aes(x=day, y= PD.mean, color = System, fill = System))+
  geom_errorbar(aes(ymin=PD.mean-PD.std, ymax=PD.mean+PD.std, group = System), size = 0.5,width=0,
                position=position_dodge(0.5), alpha = 0.6)+
  geom_point(data = alpha.PD, aes(x = day, y = PD, group = System),
             alpha = 0.6, position = position_dodge(width = 0.5), size = 2, stroke = 0)+
  geom_point(aes(x = day, y = PD.mean, group = System),
             alpha = 1, position = position_dodge(width = 0.5), size = 3)+
  theme_bw(base_size = 10) +
  scale_color_manual(values=system_color2, labels = c(expression(italic("P. inhibens (WT)")),expression(Delta*"tdaB"), "Control"))+
  scale_fill_manual(values=system_color2, labels = c(expression(italic("P. inhibens (WT)")),expression(Delta*"tdaB"), "Control"))+
  labs(y = "Faith's PD\n", x="\nTime (days)")+
  theme(strip.text.x = element_blank()) +
  theme(legend.text = element_text(colour = "black"),
        legend.title = element_text(face = "bold",colour = "black"),
        legend.position= c(0.9,0.9))+
  theme(axis.text = element_text(color = "Black", face = "bold"),
        axis.title = element_text(color = "Black", face = "bold"))+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
Fig3B 


Fig3AB = ggarrange(Fig3A,Fig3B, 
          labels = c("C", "D"),
          ncol = 2, nrow = 1, common.legend = TRUE)

ggsave("Fig3AB.tiff",
       plot = Fig3AB,
       units="cm",
       width=20,
       height=8,
       dpi = 400)


# Fig3AB = ggarrange(Fig3A + rremove("x.text"),Fig3B, 
#                    labels = c("C", "D"),
#                    ncol = 2, nrow = 1, common.legend = TRUE)


#############
### Stats ###
#############

# Linear mixed model (initial model with interaction)
lmm.alpha0 <- lmer(PD ~ System *day +(1|lineage), data=alpha.PD)

# Check studentized residuals for outliers
hist(rstudent(lmm.alpha0))

# Find values with studentized residuals above 2.5 (above 2 is a reasonable cutoff)
highResIndx=which(abs(rstudent(lmm.alpha0))>2.5)

# Linear mixed model with interaction and removed outliers
lmm.alpha <- lmer(PD ~ System *day +(1|lineage), data=alpha.PD[-highResIndx,])

# Should random effect be included
summary(lmm.alpha) # lineage has an effect so should be a random effect (not zero)

#Test model
Anova(lmm.alpha) #interaction is borderline, time is significant

#Post-hoc analysis
emm.alpha = emmeans(lmm.alpha, ~ System| day )
contrast(emm.alpha, interaction = "pairwise", adjust = "Bonferroni") #day 42 has significance vs TDA

# Get measures 
sum.alphadiv = alphadiv %>%
  group_by(day, System) %>%
  summarise(mean.observed = mean(Observed),
            sd.observed = sd(Observed),
            count = n())

