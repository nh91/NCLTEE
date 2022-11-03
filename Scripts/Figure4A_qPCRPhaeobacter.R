#QPCR - Phaeobacter inhibens

#Data import
qpcr <- read_excel("Data/qpcr_algaestd_aug2021.xlsx")
str(qpcr)

#Classes of data 
qpcr$BioRep = as.factor(qpcr$BioRep)
qpcr$Timepoint = as.factor(qpcr$Timepoint)
qpcr$TechRep = as.factor(qpcr$TechRep)
qpcr$Sample_ID = as.character(qpcr$Sample_ID)
qpcr$ExpCondition = as.factor(qpcr$ExpCondition)
str(qpcr)

#Remove D1 and D3, due to failed DNA extraction
qpcr = qpcr[qpcr$Sample_ID != "D1-0" & qpcr$Sample_ID !="D3-0",]

# Summarise data
#TechRep
qpcr_sum = qpcr %>%
  group_by(Time, ExpCondition, BioRep) %>%
  summarise( mean_copies = mean(Copies),
             sd_copies = sd(Copies),
             tech_count = n())

#BioRep 
qpcr_sum1 = qpcr_sum %>%
  group_by(Time, ExpCondition) %>%
  summarise( mean_c = mean(mean_copies),
             sd_c = sd(mean_copies),
             tech_count = n())

#Rename factor
qpcr_sum1$ExpCondition = factor(qpcr_sum1$ExpCondition, levels= c("TDA", "NoTDA", "Control"))
qpcr_sum1 = dplyr::rename(qpcr_sum1, System = ExpCondition)


#Plot 
system_color2 = c("#9C3206", "#00A091", "black")

p = ggplot(as.data.frame(qpcr_sum1), aes( x = Time, y=mean_c, color= System, shape = System)) + 
  geom_line( size = 1)+
  geom_point( size = 2)+
  geom_ribbon(aes(ymax = mean_c + sd_c, ymin = mean_c - sd_c, fill = System), alpha = 0.2, color = NA) +
  #facet_grid(exp.condition~.)+
  labs(x="Time (days)",y=" log10(CFU/mL)") +
  theme_bw(base_size = 10) + 
  scale_fill_manual(values=system_color2, labels = c(expression(italic("P. inhibens (WT)")),expression(Delta*"tdaB"), "Control"))+
  scale_color_manual(values=system_color2, labels = c(expression(italic("P. inhibens (WT)")),expression(Delta*"tdaB"), "Control"))+
  scale_y_continuous(limits = c(0,9), breaks = c(0,1,2,3,4,5,6,7,8,9))+
  scale_x_continuous(breaks = c(0,7,14,21,28,35,42,49,56,63,70))+
  theme(axis.text.y = element_text(colour = "black", face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold"), 
        legend.text = element_text(colour ="black"),
        legend.position = "none", axis.title.y = element_text(face = "bold"), 
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title.x = element_text(face = "bold", colour = "black"), 
        legend.title = element_text(colour = "black", face = "bold"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))

p
ggsave("Fig5/Fig5A.pdf",
       plot = p,
       units="cm",
       width=7,
       height=7,
       dpi = 300)


#Get legends
legend.Fig6AB.names <- as_ggplot(ggpubr::get_legend(p))

ggsave("legend.names.fig6AB.pdf",
       plot = legend.Fig6AB.names,
       units="cm",
       width=4,
       height=4,
       dpi = 300)



# ggsave(filename = "Figures/2B_QPCR.pdf",plot = p,width=12,height=8,units="cm")
# ggsave(filename = "Figures/2B_QPCR.svg",plot = p,width=12,height=8,units="cm", device = "svg")
# ggsave(filename = "Figures/2B_QPCR.png",plot = p,width=12,height=8,units="cm", device = "png")

# zoom in on day 63-70 on TDA and NoTDA
qpcr_sum$ExpCondition = factor(qpcr_sum$ExpCondition, levels= c("TDA", "NoTDA", "Control"))
qpcr_sum = dplyr::rename(qpcr_sum, System = ExpCondition)
Fig6A.zoom = ggplot(as.data.frame(qpcr_sum1) %>% filter(49<Time) %>% filter(System != "Control"),
                   aes(x = as.factor(Time), y = mean_c, fill=System, color = System, shape = System)) +
  geom_point(data = qpcr_sum %>% filter(49<Time) %>% filter(System != "Control"),
             aes(x = as.factor(Time), y = mean_copies, fill=System, stroke = 0, shape = System),
             alpha = 0.6, stat ="identity", position=position_dodge(0.5))+
  geom_point(alpha = 1, stat ="identity", position=position_dodge(0.5), size = 2)+
  geom_errorbar(aes(x = as.factor(Time), ymax = mean_c + sd_c, ymin =  mean_c - sd_c, 
                    color = System, group = System), size = 0.5, width = 0, 
                position=position_dodge(0.5), alpha = 0.6)+
  theme_bw(base_size = 10) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        strip.text = element_blank(), 
        legend.title = element_text(face = "bold", colour ="black"),
        legend.text = element_text(face = "bold", colour ="black"), 
        legend.position = "none",
        axis.text.x = element_text(face = "bold", colour ="black"),
        axis.title = element_text(face = "bold", colour ="black"),
        axis.text = element_text(face = "bold", colour ="black")) +
  scale_fill_manual(values = c("#9C3206","#00A091"))+
  scale_color_manual(values = c("#9C3206","#00A091")) +
  scale_y_continuous(limits = c(0,5), breaks = c(0,1,2,3,4,5))+
  labs(x = "", y = "")
Fig6A.zoom

ggsave("Fig5/Fig5A.zoom.pdf",
       plot = Fig6A.zoom,
       units="cm",
       width=3.5,
       height=3.5,
       dpi = 300)

#############
### Stats ###
#############

#Make factors ready
qpcr_sum$day = as.factor(qpcr_sum$Time)
qpcr_sum$lineage = paste(qpcr_sum$ExpCondition,qpcr_sum$BioRep, sep = "")

#Outlier check
lmm.qpcr <- lmer(mean_copies ~ System*day + (1|lineage), data=qpcr_sum)
Anova(lmm.qpcr)
plot(lmm.qpcr)
qqnorm(resid(lmm.qpcr))
hist(rstudent(lmm.qpcr))
highResIndx=which(abs(rstudent(lmm.qpcr))>2.5)

#Linear mixed model
lmm.qpcr <- lmer(mean_copies ~ System*day + (1|lineage), data=qpcr_sum[-highResIndx,])
#lm = lm(mean_copies ~ ExpCondition*day, data=qpcr_sum)
summary(lmm.qpcr) # lineage does have an effect so should be a random effect 
Anova(lmm.qpcr)
emm.qpcr = emmeans(lmm.qpcr, ~ System | day)
test = contrast(emm.qpcr, interaction = "pairwise", adjust = "Bonferroni")
write.csv(test, file = "Fig5/test.qpcr.csv")

str(emm.qpcr)



###########################
### Plot with dilutions ###
###########################
test1 = as.data.frame(qpcr_sum1)
test = test1[test1$Time == '14' | test1$Time == '28'| test1$Time == '42' | test1$Time == '56',]
test$mean_c = log10((10^(test$mean_c))/10)
test$Time = test$Time+0.1


newsystem = rbind(test, test1)
newsystem$System = paste(newsystem$System,"Dilution", sep ="_")
newsystem$sd_c = 0

qpcr_dilution = rbind(qpcr_sum1, newsystem)


system_color3 = c( "#5B1A18", "#7294D4", "#D8A499", "#D3D3D3", "#D3D3D3", "#D3D3D3")
qpcr_dilution$System = factor(qpcr_dilution$System, levels= c("TDA_Dilution", "NoTDA_Dilution", "Control_Dilution", "TDA", "NoTDA", "Control"))

p2 = ggplot(qpcr_dilution, aes( x = Time, y=mean_c, color= System, shape = System)) + 
  geom_line( size = 1)+
  geom_point( size = 3)+
  geom_ribbon(aes(ymax = mean_c + sd_c, ymin = mean_c - sd_c, fill = System), alpha = 0.2, color = NA) +
  #facet_grid(exp.condition~.)+
  labs(x="\n Time (days) \n",y=" \n log10(cells/mL)") +
  theme_bw(base_size = 8) + 
  scale_color_manual(values=system_color3)+
  scale_fill_manual(values=system_color3)+
  scale_shape_manual(values=c(0, 1, 2, 15, 16, 17))+
  #scale_y_continuous(breaks = c(0,1,2,3,4,5,6))+
  #scale_x_continuous(breaks = c(0,7,14,21,28,35,42,49,56,63,70))+
  theme(axis.text.y = element_text(colour = "black", face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold"), 
        legend.text = element_text(face ="bold", colour ="black"),
        legend.position = "none", axis.title.y = element_text(face = "bold"), 
        plot.title = element_text(face = "bold", hjust = 0.5),
        axis.title.x = element_text(face = "bold", colour = "black"), 
        legend.title = element_text(colour = "black", face = "bold"), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1))

p2

