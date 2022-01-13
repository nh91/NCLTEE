##########################################################

                  ### Alphadiversity ###

#########################################################

#################
### Load data ###
#################

load("Data/ps.genus.reduced.RData")


#########################
### Observed richness ###
#########################

#Make alpha estimation (Observed)
set.seed(41)
observed = estimate_richness(ps.genus.reduced, measures = "Observed")
# Make df
# Get sample ID as rows
rownames(observed) = sample_names(ps.genus.reduced) 

# Extract sample_data fra PS
s <- data.frame(sample_data(ps.genus.reduced)) 

# Get a vector with SampleID - for merging
SampleID <- row.names(observed)

#combine the observed with the sampleID column
even_stats <- data.frame(SampleID, observed) 

# merge metadata and observed  by row.names (by = 0)
alphadiv <- merge(even_stats, s, by = 0) 

# Rename factor, order level, add day and lineage
alphadiv = alphadiv %>% rename(System = exp.condition)
alphadiv$System = factor(alphadiv$System, levels= c("TDA", "NoTDA", "Control"))
alphadiv$day = factor(alphadiv$day, levels = c("0", "7", "14", "21", "28", "35", "42", "49", "56", "63", "70"))
alphadiv$lineage = paste(alphadiv$System,alphadiv$bio.rep, sep = "")

# Summarize
alphadiv_stats = alphadiv %>%
  group_by(time, System) %>%
  summarise(count = n(),
            observed.mean = mean(Observed),
            observed.std = sd(Observed))

table.observed = alphadiv_stats %>%
  arrange(System, time) %>%
  rename(Day = time)


##############
### Figure ###
##############

# color palette
system_color2 = c("#5B1A18", "#7294D4", "#D8A499")

# Figure
Fig2 = ggplot(alphadiv, aes(x=as.factor(day), y= Observed, color = System, fill = System))+
  geom_boxplot(outlier.colour = NA, color = "black", alpha = 0.9, position=position_dodge())+
  facet_grid(.~day, scale = "free_x")+
  geom_jitter(position=position_dodge(0.75), color ="black", shape = 21, aes(fill = System))+
  theme_bw(base_size = 8) +
  scale_y_continuous(limits=c(15,35), breaks=seq(15,35,by =5))+
  scale_color_manual(values=system_color2)+
  scale_fill_manual(values=system_color2)+
  labs(y = "Observed richness\n", x="\nTime (day)")+
  theme(strip.text.x = element_blank()) +
  theme(legend.text = element_text(face = "bold",colour = "black"),
        legend.title = element_text(face = "bold",colour = "black"),
        legend.position="top")+
  theme(axis.text = element_text(color = "Black", face = "bold"),
        axis.title = element_text(color = "Black", face = "bold"))+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
Fig2

ggsave("Fig2.svg",
       plot = Fig2,
       units="cm",
       width=18,
       height=8,
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

#write.csv(sum.alphadiv, file = "observedrichness_table.csv")
