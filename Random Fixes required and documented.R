colnames(Condition_Assignment) <- c("Metadata_Row", "Metadata_Column", "Metadata_ECMs", "Metadata_Ratio")

write.csv(Condition_Assignment, "H:/Lab_members/Hannah/BA e2/BA e2/R Analysis BAe2/Conditions_Assignment.csv")

for(i in 1:nrow(ds_e_phenonocont)){
  if (ds_e_phenonocont[i, "Metadata_Ratio"] == "\"50:50\""){
    
    ds_e_phenonocont[i, "Metadata_Ratio"] <- "50:50" 
  
    
  }else if (ds_e_phenonocont[i, "Metadata_Ratio"] == "\"75:25\""){
    
    ds_e_phenonocont[i, "Metadata_Ratio"] <- "75:25"   
    
  }}



interaction.plot(ECMsPairwiseSplitUpdated$ECM_Combo_2, ECMsPairwiseSplitUpdated$ECM_Combo_1, ECMsPairwiseSplitUpdated$p.value, las=2, col = 1:29, trace.label = "ECM Combo 1", ylab = "p.value", xlab = " ")

ECMSPairwisels = as.data.frame(ECMSPairwise$lsmeans)

plot(ECMSPairwisels$Metadata_ECMs, ECMSPairwisels$lsmean)

colors <- c(ECMSPairwisels$Metadata_ECMs) 

ggplot(data = ECMSPairwisels,
       aes(x = Metadata_ECMs,
           y = lsmean)) + 
  theme(axis.title.x = element_text(size= 30),
        axis.text.x = element_text(angle = 90, size= 15),
        axis.title.y = element_text(size= 30),
        axis.text.y = element_text(size= 15),
        strip.text =  element_text(size = 30),
        plot.title = element_text(size = 30))+
  scale_fill_manual()
  geom_col() +
  xlab("LS Mean") +
  ylab("ECMs")
  
require(ggpubr)
require(ggplot2)  
  
ggplot(data = ECMSPairwisels)  + 
  geom_segment( aes(x = Metadata_ECMs, xend = Metadata_ECMs, y = lower.CL, yend = upper.CL), color = "blue") +
  geom_point( aes(x = Metadata_ECMs, y = lower.CL), color="darkblue", shape = 1) + 
  geom_point( aes(x = Metadata_ECMs, y = upper.CL), color="darkblue", shape = 1) + 
  geom_point( aes(x = Metadata_ECMs, y = lsmean), color="black", shape = 4) + 
  xlab("ECM Conditions")+ 
  ylab("Least Square Mean") + 
  #stat_pvalue_manual(ECMsPvaluestatsig, label = "p.value", y.position = 6, step.increase = 0.1, hide.ns = TRUE) +
  coord_flip()
  

ECMsPvaluestat <- subset(ECMsPairwiseSplitUpdated, select = c("ECM_Combo_1", "ECM_Combo_2", "p.value"))
colnames(ECMsPvaluestat) <- c("group1", "group2", "p.value")
ECMsPvaluestatsig <- subset(ECMsPvaluestat, p.value <= 0.05)

ggplot(data = ECMsPvaluestat,
        aes(x = group1, y = group2, fill = p.value))+
  geom_tile() +
  theme(axis.title.x = element_text(size= 30),
        axis.text.x = element_text(angle = 90, size= 15),
        axis.title.y = element_text(size= 30),
        axis.text.y = element_text(size= 15),
        strip.text =  element_text(size = 30),
        plot.title = element_text(size = 30))


## ECM:Ratio
require(ggpubr)
require(ggplot2)  

ECMSRatioPairwisels = as.data.frame(ECMbyRatioPair$lsmeans)

ggplot(data = ECMSRatioPairwisels,
       aes(x = Metadata_ECMs,
           y = lsmean)) + 
  theme(axis.title.x = element_text(size= 20),
        axis.text.x = element_text(angle = 0, size= 15),
        axis.title.y = element_text(size= 20),
        axis.text.y = element_text(size= 10),
        strip.text =  element_text(size = 20),
        plot.title = element_text(size = 30))+
  scale_fill_manual()+
  facet_wrap("Metadata_Ratio") +
geom_col() +
  xlab("LS Mean") +
  ylab("ECMs")+
  coord_flip()

ECMRatiols25 <- subset(ECMSRatioPairwisels, ECMSRatioPairwisels$Metadata_Ratio == "25:75")
ECMRatiols50 <- subset(ECMSRatioPairwisels, ECMSRatioPairwisels$Metadata_Ratio == "50:50")
ECMRatiols75 <- subset(ECMSRatioPairwisels, ECMSRatioPairwisels$Metadata_Ratio == "75:25")

ggplot(data = ECMRatiols25,
       aes(x = Metadata_ECMs,
           y = lsmean)) + 
  theme(axis.title.x = element_text(size= 30),
        axis.text.x = element_text(angle = 90, size= 15),
        axis.title.y = element_text(size= 30),
        axis.text.y = element_text(size= 15),
        strip.text =  element_text(size = 30),
        plot.title = element_text(size = 30))+
  scale_fill_manual()+
  #facet_wrap("Metadata_Ratio") +
  geom_col() +
  xlab("LS Mean") +
  ylab("ECMs") + 
  title("25:75")+
  coord_flip()

ggplot(data = ECMSRatioPairwisels)  + 
  geom_segment( aes(x = Metadata_ECMs, xend = Metadata_ECMs, y = lower.CL, yend = upper.CL), color = "blue") +
  geom_point( aes(x = Metadata_ECMs, y = lower.CL), color="darkblue", shape = 1) + 
  geom_point( aes(x = Metadata_ECMs, y = upper.CL), color="darkblue", shape = 1) + 
  geom_point( aes(x = Metadata_ECMs, y = lsmean), color="black", shape = 4) + 
  facet_wrap("Metadata_Ratio")+
  xlab("ECM Conditions")+ 
  ylab("Least Square Mean") + 
  #stat_pvalue_manual(ECMsPvaluestatsig, label = "p.value", y.position = 6, step.increase = 0.1, hide.ns = TRUE) +
  coord_flip()  

ECMsRatioPvaluestat <- subset(ECMRatioSplit, select = c("ECM_Combo_1", "ECM_Combo_2", "Metadata_Ratio", "p.value"))
ECMsRatioPvalue25 <- subset(ECMsRatioPvaluestat, ECMsRatioPvaluestat$Ratio == "25:75")
colnames(ECMsRatioPvaluestat) <- c("group1", "group2", "Ratio", "p.value")
ECMsRatioPvaluestatsig <- subset(ECMsRatioPvaluestat, p.value <= 0.05)

ggplot(data = ECMsRatioPvalue25,
       aes(x = group1, y = group2, fill = p.value))+
  geom_tile() 
  # theme(axis.title.x = element_text(size= 30),
  #       axis.text.x = element_text(angle = 90, size= 15),
  #       axis.title.y = element_text(size= 30),
  #       axis.text.y = element_text(size= 15),
  #       strip.text =  element_text(size = 30),
  #       plot.title = element_text(size = 30))

for(i in 1:nrow(X4_21_22HM)){
  if (X4_21_22HM[i, "Metadata_Ratio"] == "25:75"){
    
    X4_21_22HM[i, "Metadata_Ratio"] <- "25:75"
    
  }else if (X4_21_22HM[i, "Metadata_Ratio"] == "\"50:50\""){
    
    X4_21_22HM[i, "Metadata_Ratio"] <- "50:50"  
    
  }else if (X4_21_22HM[i, "Metadata_Ratio"] == "\"75:25\""){
    
    X4_21_22HM[i, "Metadata_Ratio"] <- "75:25" }}

Day2nohundo$Metadata_Ratio <- as.factor(Day2nohundo$Metadata_Ratio)
Day2nohundo$Metadata_Ratio <- factor(Day2nohundo$Metadata_Ratio, levels=c("25:75", "\"50:50\"", "\"75:25\""))

HM4_21_22 <- ggplot(data = X4_21_22HM,
                                aes(x = Metadata_Ratio,
                                    y = Metadata_ECMs_Abb,
                                    fill = Mean_Object_Counts)) +
  facet_wrap(vars(Metadata_Stiffness)) +
  theme(axis.title.x = element_text(size= 40),
        axis.text.x = element_text(angle = 0, size= 20),
        axis.title.y = element_text(size= 40),
        axis.text.y = element_text(size= 20),
        strip.text =  element_text(size = 30),
        plot.title = element_text(size = 30),
        legend.text=element_text(size=20),
        legend.title = element_text(size=20))+
  geom_tile() +
  xlab("Ratio") +
  ylab("ECMs")+
  labs(title = "Number of Cells on Islands Day 2", fill = "Number of Cells")
ggsave(filename = "plots/4_21_22HM2.png",
       width = 10, height= 10)

