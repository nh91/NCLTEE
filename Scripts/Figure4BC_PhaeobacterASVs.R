######################################################################################

              ### Relative abundance of Phaeobacter ASVs ###

######################################################################################


#################
### Load data ###
#################

#PS object with ASVs
load("Data/ReadO/ps.asv.reduced.wTree.RData")

########################
### Phaeobacter ASVs ###
########################

# PS object with Phaeobacter ASVs
ps.phaeobacter = ps.new %>% 
  phyloseq::transform_sample_counts(function(x) {x/sum(x)}*100 ) %>%
  phyloseq::subset_taxa(genus == "Phaeobacter") %>%
  phyloseq::prune_taxa(taxa_sums(.) >0, .)

# Get reference sequences for BLAST 
refseq.phaeobacter = as.data.frame(refseq(ps.phaeobacter))
#write.csv(refseq.phaeobacter, file = "Data/refseq.phaeobacter.csv")


######################
### Data wrangling ###
######################

# Make df for plotting
df.phaeobacter = ps.phaeobacter %>%
  psmelt() %>%
  dplyr::rename(System = exp.condition)
df.phaeobacter$OTU = as.factor(df.phaeobacter$OTU)

write_csv(df.phaeobacter, file = "ASVabundances.csv")
#Calculate mean for bio.rep
sum.phaeobacter = df.phaeobacter %>%
  group_by(day, System, OTU) %>%
  summarise(mean.abundance = mean(Abundance),
            sd.abundance = sd(Abundance),
            count = n())

##############
### Figure ###
##############

# Change levels
sum.phaeobacter$System = factor(sum.phaeobacter$System , levels= c("TDA", "NoTDA", "Control"))
system_color2 = c("#9C3206", "#00A091", "black")
labels<- c(expression(italic("P. inhibens (WT)")),expression(Delta*"tdaB"), "Control")
names(labels) = c("TDA", "NoTDA", "Control")

# Line plot
Fig6C = 
  ggplot(sum.phaeobacter, 
               aes(x = day, y = mean.abundance, 
                   fill=OTU, color = OTU, shape = System)) +
  geom_line(size = 1)+
  geom_point(size = 2)+
  geom_ribbon(aes(ymax = mean.abundance + sd.abundance, 
                  ymin = mean.abundance - sd.abundance, fill = OTU), alpha = 0.2, color = NA) +
  facet_wrap(vars(System), nrow = 3) + 
  theme_bw(base_size = 10) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        strip.background = element_rect(fill = "white"),
        legend.position = "none",
        strip.text = element_text(face = "bold", colour ="black"), 
        legend.title = element_text(face = "bold", colour ="black"),
        legend.text = element_text(face = "bold", colour ="black"),
        axis.text.x = element_text(face = "bold", colour ="black"),
        axis.title = element_text(face = "bold", colour ="black"),
        axis.text = element_text(face = "bold", colour ="black")) +
  scale_fill_manual(values = c("#CE738C", "Black"))+
  scale_color_manual(values = c( "#CE738C", "Black"))+
  labs(x = "Time (days)", y = " Relative abundance (%)") + 
  scale_y_continuous(breaks = c(0,2,4,6,8,10,12,14))+ 
  scale_x_continuous(breaks = c(0,7,14,21,28,35,42,49,56,63,70))
Fig6C

#Save figure
ggsave("Fig5/Fig5B.pdf",
       plot = Fig6C,
       units="cm",
       width=8,
       height=18,
       dpi = 300)

legend.Fig6C <- cowplot::get_legend(Fig6C)
grid.newpage()
grid.draw(legend.Fig6C)

# zoom in on day 42-70
df.phaeobacter$System = factor(df.phaeobacter$System , levels= c("TDA", "NoTDA", "Control"))
Fig6C.con = ggplot(sum.phaeobacter %>% filter(35<day) %>% filter(System == "Control"),
                aes(x = day, y = mean.abundance, fill=OTU, color = OTU)) +
  geom_point(alpha = 1, stat ="identity", position=position_dodge(3), size = 2, shape = 3)+
  geom_point(data = df.phaeobacter %>% filter(35<day) %>% filter(System == "Control"),
             aes(x = day, y = Abundance,fill=OTU),shape = 15, size = 2,
               alpha = 0.8, stat ="identity", position=position_dodge(3))+
  geom_errorbar(aes(x = day, ymax = mean.abundance + sd.abundance, ymin = mean.abundance - sd.abundance, 
                    color = OTU, group = OTU), size = 0.5, width = 0, position=position_dodge(3), alpha = 0.9)+
  theme_bw(base_size = 10) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text = element_blank(), 
        legend.title = element_text(face = "bold", colour ="black"),
        legend.text = element_text(face = "bold", colour ="black"), 
        legend.position = "",
        axis.text.x = element_text(face = "bold", colour ="black"),
        axis.title = element_text(face = "bold", colour ="black"),
        axis.text = element_text(face = "bold", colour ="black")) +
  scale_fill_manual(values = c("#CE738C", "Black"))+
  scale_color_manual(values = c( "#CE738C", "Black")) +
  labs(x = "", y = "") + 
  scale_y_continuous(limits = c(-0.1,2), breaks = c(0,0.5,1,1.5, 2.0, 2.5))+
  scale_x_continuous(breaks = c(42,49,56,63,70))
Fig6C.con

ggsave("Fig5/fig5B.con.pdf",
       plot = Fig6C.con,
       units="cm",
       width=5,
       height=5,
       dpi = 300)

#############
### Stats ###
#############

#Make df for TDA
df.TDA = df.phaeobacter[df.phaeobacter$System == "TDA",] %>% select(Abundance, System, bio.rep, day, OTU)
df.TDA$lineage = paste(df.TDA$System,df.TDA$bio.rep, sep = "")
df.TDA$day = as.factor(df.TDA$day)

#Make df for NoTDA
df.nTDA = df.phaeobacter[df.phaeobacter$System == "NoTDA",] %>% select(Abundance, System, bio.rep, day, OTU)
df.nTDA$lineage = paste(df.nTDA$System,df.nTDA$bio.rep, sep = "")
df.nTDA$day = as.factor(df.nTDA$day)

#Make df for Control
df.con = df.phaeobacter[df.phaeobacter$System == "Control",] %>% select(Abundance, System, bio.rep, day, OTU)
df.con$lineage = paste(df.con$System,df.con$bio.rep, sep = "")
df.con$day = as.factor(df.con$day)

#TDA
#Linear mixed model for ASV difference in abundance in TDA
lmm.alpha0 <- lmer(Abundance ~ OTU*day +( 1|lineage), data=df.TDA)
qqnorm(resid(lmm.alpha0))
plot(lmm.alpha0)
# Check studentized residuals for outliers
hist(rstudent(lmm.alpha0))
# Find values with studentized residuals above 2.5 (above 2 is a reasonable cutoff)
highResIndx=which(abs(rstudent(lmm.alpha0))>2.5)

#Model without outliers
lmm.TDA <- lmer(Abundance ~ OTU*day +( 1|lineage), data=df.TDA[-highResIndx,])
qqnorm(resid(lmm.TDA))
plot(lmm.TDA)
summary(lmm.TDA) # lineage does have an effect so should be a random effect 
Anova(lmm.TDA)
emm.TDA = emmeans(lmm.TDA, ~ OTU | day)
lmmTDA = data.frame(contrast(emm.TDA, interaction = "pairwise", adjust = "Bonferroni"))
lmmTDA$System = "P. inhibens (WT)"

#Linear mixed model for ASV difference in abundance in nTDA
lmm.alpha0 <- lmer(Abundance ~ OTU*day +( 1|lineage), data=df.nTDA)
qqnorm(resid(lmm.alpha0))
plot(lmm.alpha0)
# Check studentized residuals for outliers
hist(rstudent(lmm.alpha0))
# Find values with studentized residuals above 2.5 (above 2 is a reasonable cutoff)
highResIndx=which(abs(rstudent(lmm.alpha0))>2.5)

#Model without outliers
lmm.nTDA <- lmer(Abundance ~ OTU*day +( 1|lineage), data=df.nTDA[-highResIndx,])
qqnorm(resid(lmm.nTDA))
plot(lmm.nTDA)
summary(lmm.nTDA) # lineage does have an effect so should be a random effect 
Anova(lmm.nTDA)
emm.nTDA = emmeans(lmm.nTDA, ~ OTU | day)
lmmnTDA = data.frame(contrast(emm.nTDA, interaction = "pairwise", adjust = "Bonferroni"))
lmmnTDA$System = "delta-tdaB"

#Linear mixed model for ASV difference in abundance in Control
lmm.alpha0 <- lmer(Abundance ~ OTU*day +( 1|lineage), data=df.con)
qqnorm(resid(lmm.alpha0))
plot(lmm.alpha0)
# Check studentized residuals for outliers
hist(rstudent(lmm.alpha0))
# Find values with studentized residuals above 2.5 (above 2 is a reasonable cutoff)
highResIndx=which(abs(rstudent(lmm.alpha0))>2.5)
lmm.con <- lmer(Abundance ~ OTU*day +( 1|lineage), data=df.con[-highResIndx,])
summary(lmm.con) # lineage does have an effect so should be a random effect 
Anova(lmm.con)
emm.con = emmeans(lmm.con, ~ OTU | day)
lmmcon = data.frame(contrast(emm.con, interaction = "pairwise", adjust = "Bonferroni"))
lmmcon$System = "Control"

lmmASVSystem = rbind(lmmTDA,lmmnTDA, lmmcon)
#write.csv(lmmASVSystem, file = "Fig5/lmmASVSystem.csv")



# Get 42 to 70 for each asv
Fig6D.ASVs = ggplot(sum.phaeobacter %>% filter(35<day),
                   aes(x = day, y = mean.abundance, fill=System, color = System)) +
  geom_point(data = df.phaeobacter %>% filter(35<day),
             aes(x = day, y = Abundance,fill=System, stroke = 0, shape = System),
             alpha = 0.6, stat ="identity", position=position_dodge(3))+
  geom_point(aes(shape = System),alpha = 1, stat ="identity", position=position_dodge(3), size = 2.5)+
  geom_errorbar(aes(x = day, ymax = mean.abundance + sd.abundance, ymin = mean.abundance - sd.abundance, 
                    color = System, group = System), size = 0.5, width = 0, position=position_dodge(3), alpha = 0.6)+
  facet_wrap(vars(OTU), nrow = 2)+
  theme_bw(base_size = 10) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        legend.title = element_text(face = "bold", colour ="black"),
        legend.text = element_text(face = "bold", colour ="black"), 
        strip.background = element_rect(fill = "white"),
        legend.position = "none",
        strip.text = element_text(face = "bold", colour ="black"), 
        axis.text.x = element_text(face = "bold", colour ="black"),
        axis.title = element_text(face = "bold", colour ="black"),
        axis.text = element_text(face = "bold", colour ="black")) +
  scale_fill_manual(values = system_color2 )+
  scale_color_manual(values = system_color2 ) +
   labs(x = "Time (days)", y = "Relative abundance (%)") + 
  scale_y_continuous(limits = c(-0.1,2.7), breaks = c(0,0.5,1,1.5,2.0,2.5))+
  scale_x_continuous(breaks = c(42,49,56,63,70))
Fig6D.ASVs

ggsave("Fig5/Fig5C.ASVs.pdf",
       plot = Fig6D.ASVs,
       units="cm",
       width=7,
       height=10,
       dpi = 300)

#############
### Stats ###
#############

#Make df for ASV18
df.18 = df.phaeobacter[df.phaeobacter$OTU == "ASV_18",] %>% select(Abundance, System, bio.rep, day, OTU)
df.18$lineage = paste(df.18$System,df.18$bio.rep, sep = "")
df.18$day = as.factor(df.18$day)

#Make df for ASV30
df.30 = df.phaeobacter[df.phaeobacter$OTU == "ASV_30",] %>% select(Abundance, System, bio.rep, day, OTU)
df.30$lineage = paste(df.30$System,df.30$bio.rep, sep = "")
df.30$day = as.factor(df.30$day)


#Linear mixed model for ASVS18- difference in abundance between system 
# Linear mixed model (initial model with interaction)
lmm.alpha0 <- lmer(Abundance ~ System*day +(1|lineage), data=df.18)
qqnorm(resid(lmm.alpha0))
plot(lmm.alpha0)
# Check studentized residuals for outliers
hist(rstudent(lmm.alpha0))
# Find values with studentized residuals above 2.5 (above 2 is a reasonable cutoff)
highResIndx=which(abs(rstudent(lmm.alpha0))>2.5)

lmm.18 <- lmer(Abundance ~ System*day +( 1|lineage), data=df.18[-highResIndx,])
qqnorm(resid(lmm.18))
plot(lmm.18)
summary(lmm.18) # lineage does have an effect so should be a random effect 
Anova(lmm.18)
emm.18 = emmeans(lmm.18, ~ System | day)
lmmASV18 = data.frame(contrast(emm.18, interaction = "pairwise", adjust = "Bonferroni"))
lmmASV18$ASV = "ASV18"


#Linear mixed model ASV330- difference in abundance between system 
# Linear mixed model (initial model with interaction)
lmm.alpha0 <- lmer(Abundance ~ System*day +(1|lineage), data=df.30)
qqnorm(resid(lmm.alpha0))
plot(lmm.alpha0)
# Check studentized residuals for outliers
hist(rstudent(lmm.alpha0))
# Find values with studentized residuals above 2.5 (above 2 is a reasonable cutoff)
highResIndx=which(abs(rstudent(lmm.alpha0))>2.5)

lmm.30 <- lmer(Abundance ~ System*day +( 1|lineage), data=df.30[-highResIndx,])
qqnorm(resid(lmm.30))
plot(lmm.30)
summary(lmm.30) # lineage does have an effect so should be a random effect 
Anova(lmm.30)
emm.30 = emmeans(lmm.30, ~ System | day)
contrast(emm.30, interaction = "pairwise", adjust = "Bonferroni")
lmmASV30 = data.frame(contrast(emm.30, interaction = "pairwise", adjust = "Bonferroni"))
lmmASV30$ASV = "ASV30"

lmmASVPhaeo = rbind(lmmASV18, lmmASV30)
write.csv(lmmASVPhaeo, file = "Fig5/lmmASVPhaeo.csv")
