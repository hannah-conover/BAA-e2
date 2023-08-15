require(ggplot2)
require(ggthemes)
require(plyr)
require(stargazer)
require(RColorBrewer)
require(scales)
require(gtable)
require(gridExtra)
require(data.table)
require(dplyr)
require(ggplot2)
require(readr)
require(tidyr)

# Treatment Assignment ----------------------------------------------------
blue_objects$Metadata_Batch <- "None"
blue_objects$Metadata_Time <- "None"
blue_objects$Metadata_Stiffness <- "None"

#Assign Slide-level Treatment Conditions

for(i in 1:nrow(blue_objects)){
  if (blue_objects[i, "Metadata_Slide"] == "01"){
    
    blue_objects[i, "Metadata_Stiffness"] <- "1 kPa"
    blue_objects[i, "Metadata_Batch"] <- "2"
    blue_objects[i, "Metadata_Time"] <- "Day 0"
  
  }else if (blue_objects[i,"Metadata_Slide"]== "02"){
   
    blue_objects[i, "Metadata_Stiffness"] <- "25 kPa"
    blue_objects[i, "Metadata_Batch"] <- "2"
    blue_objects[i, "Metadata_Time"] <- "Day 0"  
    
  }else if (blue_objects[i,"Metadata_Slide"]== "03"){
    
    blue_objects[i, "Metadata_Stiffness"] <- "1 kPa"
    blue_objects[i, "Metadata_Batch"] <- "2"
    blue_objects[i, "Metadata_Time"] <- "Day 0" 
    
  }else if (blue_objects[i,"Metadata_Slide"]== "04"){
    
    blue_objects[i, "Metadata_Stiffness"] <- "25 kPa"
    blue_objects[i, "Metadata_Batch"] <- "1"
    blue_objects[i, "Metadata_Time"] <- "Day 0" 
    
  }else if (blue_objects[i,"Metadata_Slide"]== "05"){
    
    blue_objects[i, "Metadata_Stiffness"] <- "25 kPa"
    blue_objects[i, "Metadata_Batch"] <- "1"
    blue_objects[i, "Metadata_Time"] <- "Day 0" 
    
  }else if (blue_objects[i,"Metadata_Slide"]== "06"){
    
    blue_objects[i, "Metadata_Stiffness"] <- "1 kPa"
    blue_objects[i, "Metadata_Batch"] <- "1"
    blue_objects[i, "Metadata_Time"] <- "Day 0" 
    
  }else if (blue_objects[i,"Metadata_Slide"]== "07"){
    
    blue_objects[i, "Metadata_Stiffness"] <- "25 kPa"
    blue_objects[i, "Metadata_Batch"] <- "2"
    blue_objects[i, "Metadata_Time"] <- "Day 2" 
    
  }else if (blue_objects[i,"Metadata_Slide"]== "08"){
    
    blue_objects[i, "Metadata_Stiffness"] <- "25 kPa"
    blue_objects[i, "Metadata_Batch"] <- "2"
    blue_objects[i, "Metadata_Time"] <- "Day 2" 
    
  }else if (blue_objects[i,"Metadata_Slide"]== "09"){
    
    blue_objects[i, "Metadata_Stiffness"] <- "1 kPa"
    blue_objects[i, "Metadata_Batch"] <- "2"
    blue_objects[i, "Metadata_Time"] <- "Day 2" 
    
  }else if (blue_objects[i,"Metadata_Slide"]== "10"){
    
    blue_objects[i, "Metadata_Stiffness"] <- "1 kPa"
    blue_objects[i, "Metadata_Batch"] <- "1"
    blue_objects[i, "Metadata_Time"] <- "Day 2" 
    
  }else if (blue_objects[i,"Metadata_Slide"]== "11"){
    
    blue_objects[i, "Metadata_Stiffness"] <- "25 kPa"
    blue_objects[i, "Metadata_Batch"] <- "1"
    blue_objects[i, "Metadata_Time"] <- "Day 2" 
    
  }else if (blue_objects[i,"Metadata_Slide"]== "12"){
    
    blue_objects[i, "Metadata_Stiffness"] <- "1 kPa"
    blue_objects[i, "Metadata_Batch"] <- "1"
    blue_objects[i, "Metadata_Time"] <- "Day 2" 
  }
}

# Check for unique values in blue_objects
unique(blue_objects$Metadata_Batch)
unique(blue_objects$Metadata_Stiffness)
unique(blue_objects$Metadata_Time)



# Merge Island Specific Conditions ----------------------------------------

blue_object = merge(blue_objects, Condition_Assignment, by = c("Metadata_Column", "Metadata_Row"))


# Summarizing the Data ----------------------------------------------------

ds_s_pheno <- ddply(blue_object,
                    .(Metadata_ECMs,
                      Metadata_Ratio,
                      Metadata_Stiffness,
                      Metadata_Batch,
                      Metadata_Slide,
                      Metadata_Experiment,
                      Metadata_Row,
                      Metadata_Column,
                      Metadata_Time,
                      Metadata_Scene
                      ),
                    summarize,
                    Blue = mean(Intensity_MeanIntensity_blue,
                                na.rm = TRUE),
                    Object_Counts = length(ObjectNumber)
)

ds_e_pheno <- ddply(ds_s_pheno,
                   .(Metadata_ECMs,
                     Metadata_Ratio,
                     Metadata_Stiffness,
                     Metadata_Slide,
                     Metadata_Batch,
                     Metadata_Experiment,
                     Metadata_Row,
                     Metadata_Column,
                     Metadata_Time
                   ),
                   summarize,
                   Mean_Blue = mean(Blue,
                               na.rm = TRUE),
                   Mean_Object_Counts = mean(Object_Counts,
                                             na.rm = TRUE),
                   SD_Object_Counts = sd(Object_Counts,
                                         na.rm = TRUE))

ds_pheno <- ddply(ds_e_pheno,
                  .(Metadata_ECMs,
                    Metadata_Ratio,
                    Metadata_Stiffness,
                    Metadata_Experiment,
                    Metadata_Batch,
                    Metadata_Row,
                    Metadata_Column,
                    Metadata_Time
                  ),
                  summarize,
                  Mean_Mean_Blue = mean(Mean_Blue,
                                   na.rm = TRUE),
                  SD_Mean_Blue = sd(Mean_Blue,
                                    na.rm = TRUE)/sqrt(NROW(Mean_Blue)),
                  Mean_Mean_Object_Counts = mean(Mean_Object_Counts,
                                            na.rm = TRUE),
                  SD_Mean_Object_Counts = sd(Mean_Object_Counts,
                                        na.rm = TRUE)/sqrt(NROW(Mean_Object_Counts)))

# BoxCox and Models (ds_s_pheno individual islands) -------------------------------------------------------
library(daewr)
library(GAD)
library(lme4)
library(MASS)

Stiffness <- as.fixed(ds_s_pheno$Metadata_Stiffness)
Batch <- as.random(ds_s_pheno$Metadata_Batch)
Ratio <- as.fixed(ds_s_pheno$Metadata_Ratio)
ECMs <- as.fixed(ds_s_pheno$Metadata_ECMs)
Time <- as.fixed(ds_s_pheno$Metadata_Time)
Row <- as.random(ds_s_pheno$Metadata_Row)
Column <- as.random(ds_s_pheno$Metadata_Column)

model <- aov(Object_Counts ~ Stiffness + Stiffness%in%Batch + Ratio +ECMs + Time + Time%in%Batch + 
               Ratio%in%Row + Ratio%in%Column + ECMs%in%Row + ECMs%in%Column+ 
               Stiffness*Ratio + Stiffness*ECMs + Stiffness*Time + Time*ECMs + Ratio*Time + Ratio*ECMs + 
               Stiffness*Ratio*Time + Stiffness*Ratio*ECMs + Stiffness*ECMs*Time + Ratio*Time*ECMs +
             Stiffness*Ratio*Time*ECMs, data = ds_s_pheno)
summary(model)
boxcox(model)

#Typically the boxcox will show a transformation of near 0 which is the log

logmodel <- aov(log(Object_Counts)~Stiffness + Stiffness%in%Batch + Ratio +ECMs + Time + Time%in%Batch + 
                  Ratio%in%Row + Ratio%in%Column + ECMs%in%Row + ECMs%in%Column+ 
                  Stiffness*Ratio + Stiffness*ECMs + Stiffness*Time + Time*ECMs + Ratio*Time + Ratio*ECMs + 
                  Stiffness*Ratio*Time + Stiffness*Ratio*ECMs + Stiffness*ECMs*Time + Ratio*Time*ECMs +
                  Stiffness*Ratio*Time*ECMs, data = ds_s_pheno)

summary(logmodel)
anova(logmodel)
gad(logmodel)

logmodelcombo <- aov(log(Object_Counts)~Stiffness + Stiffness*Time%in%Batch + Ratio +ECMs + Time + Time%in%Batch + 
                       Ratio%in%Row + Ratio%in%Column + ECMs%in%Row + ECMs%in%Column+ 
                       Stiffness*Ratio + Stiffness*ECMs + Stiffness*Time + Time*ECMs + Ratio*Time + Ratio*ECMs + 
                       Stiffness*Ratio*Time + Stiffness*Ratio*ECMs + Stiffness*ECMs*Time + Ratio*Time*ECMs +
                       Stiffness*Ratio*Time*ECMs, data = ds_s_pheno)

wplogmodel <- aov(log(Object_Counts)~ Stiffness + Time + Stiffness*Time + Batch, data=ds_s_pheno)
summary(wplogmodel)

#Pull out Control (100% Ratios)
controls <- subset(ds_s_pheno, Metadata_Ratio == 100, select = c(Metadata_Batch, Metadata_Stiffness, Metadata_Time, Metadata_Slide, Metadata_Scene, Metadata_ECMs, Object_Counts, Metadata_Row, Metadata_Column))
ds_s_phenonocont <- subset(ds_s_pheno, Metadata_Ratio != 100, select = c(Metadata_Ratio, Metadata_Batch, Metadata_Stiffness, Metadata_Time, Metadata_Slide, Metadata_Scene, Metadata_ECMs, Object_Counts, Metadata_Row, Metadata_Column))

StiffnessC <- as.fixed(controls$Metadata_Stiffness)
BatchC <- as.random(controls$Metadata_Batch)
ECMsC <- as.fixed(controls$Metadata_ECMs)
TimeC <- as.fixed(controls$Metadata_Time)
RowC <- as.random(controls$Metadata_Row)
ColumnC <- as.random(controls$Metadata_Column)

contlogmod <- aov(log(Object_Counts)~StiffnessC + StiffnessC%in%BatchC +ECMsC + TimeC + TimeC%in%BatchC + StiffnessC*TimeC%in%BatchC + 
                     ECMsC%in%RowC + ECMsC%in%ColumnC+ 
                     StiffnessC*ECMsC + StiffnessC*TimeC + TimeC*ECMsC + 
                     StiffnessC*ECMsC*TimeC, data = controls)
summary(contlogmod)

StiffnessNH <- as.fixed(ds_s_phenonocont$Metadata_Stiffness)
BatchNH <- as.random(ds_s_phenonocont$Metadata_Batch)
SceneNH <- as.random(ds_s_phenonocont$Metadata_Scene)
RatioNH <- as.fixed(ds_s_phenonocont$Metadata_Ratio)
ECMsNH <- as.fixed(ds_s_phenonocont$Metadata_ECMs)
TimeNH <- as.fixed(ds_s_phenonocont$Metadata_Time)
RowNH <- as.random(ds_s_phenonocont$Metadata_Row)
ColumnNH <- as.random(ds_s_phenonocont$Metadata_Column)

no100logmod <- aov(log(Object_Counts)~StiffnessNH + StiffnessNH%in%BatchNH + RatioNH +ECMsNH + TimeNH + TimeNH%in%BatchNH + TimeNH*StiffnessNH%in%BatchNH +
                     RatioNH%in%RowNH + RatioNH%in%ColumnNH + ECMsNH%in%RowNH + ECMsNH%in%ColumnNH+ SceneNH +
                     StiffnessNH*RatioNH + StiffnessNH*ECMsNH + StiffnessNH*TimeNH + TimeNH*ECMsNH + RatioNH*TimeNH + RatioNH*ECMsNH + 
                     StiffnessNH*RatioNH*TimeNH + StiffnessNH*RatioNH*ECMsNH + StiffnessNH*ECMsNH*TimeNH + RatioNH*TimeNH*ECMsNH +
                     StiffnessNH*RatioNH*TimeNH*ECMsNH, data = ds_s_phenonocont)
summary(no100logmod, digits=6)


# ds_e_pheno and Box Cox --------------------------------------------------

library(daewr)
library(GAD)
library(lme4)
library(MASS)

EStiffness <- as.fixed(ds_e_pheno$Metadata_Stiffness)
EBatch <- as.random(ds_e_pheno$Metadata_Batch)
ERatio <- as.fixed(ds_e_pheno$Metadata_Ratio)
EECMs <- as.fixed(ds_e_pheno$Metadata_ECMs)
ETime <- as.fixed(ds_e_pheno$Metadata_Time)
# ERow <- as.random(ds_e_pheno$Metadata_Row)
# EColumn <- as.random(ds_e_pheno$Metadata_Column)

Emodel <- aov(Mean_Object_Counts ~ EStiffness + EStiffness%in%EBatch + ERatio +EECMs + ETime + ETime%in%EBatch + EStiffness*ETime%in%EBatch+
               # Ratio%in%Row + Ratio%in%Column + ECMs%in%Row + ECMs%in%Column+ 
               EStiffness*ERatio + EStiffness*EECMs + EStiffness*ETime + ETime*EECMs + ERatio*ETime + ERatio*EECMs + 
               EStiffness*ERatio*ETime + EStiffness*ERatio*EECMs + EStiffness*EECMs*ETime + ERatio*ETime*EECMs +
               EStiffness*ERatio*ETime*EECMs, data = ds_e_pheno)
summary(Emodel)
boxcox(Emodel)

#Typically the boxcox will show a transformation of near 0 which is the log

logEmodel <- aov(log(Mean_Object_Counts)~EStiffness + EStiffness%in%EBatch + ERatio +EECMs + ETime + ETime%in%EBatch + EStiffness*ETime%in%EBatch+
                   # Ratio%in%Row + Ratio%in%Column + ECMs%in%Row + ECMs%in%Column+ 
                   EStiffness*ERatio + EStiffness*EECMs + EStiffness*ETime + ETime*EECMs + ERatio*ETime + ERatio*EECMs + 
                   EStiffness*ERatio*ETime + EStiffness*ERatio*EECMs + EStiffness*EECMs*ETime + ERatio*ETime*EECMs +
                   EStiffness*ERatio*ETime*EECMs, data = ds_e_pheno)

summary(logEmodel)
anova(logmodel)
gad(logmodel)

wplogEmodel <- aov(log(Mean_Object_Counts)~ EStiffness + ETime + EStiffness*ETime + EBatch, data=ds_e_pheno)
summary(wplogEmodel)

#Pull out Control (100% Ratios)
Econtrols <- subset(ds_e_pheno, Metadata_Ratio == 100, select = c(Metadata_Batch, Metadata_Stiffness, Metadata_Time, Metadata_Slide, Metadata_ECMs, Mean_Object_Counts, Metadata_Row, Metadata_Column))
ds_e_phenonocont <- subset(ds_e_pheno, Metadata_Ratio != 100, select = c(Metadata_Ratio, Metadata_Batch, Metadata_Stiffness, Metadata_Time, Metadata_Slide,  Metadata_ECMs, Mean_Object_Counts, Metadata_Row, Metadata_Column))

StiffnessEC <- as.fixed(Econtrols$Metadata_Stiffness)
BatchEC <- as.random(Econtrols$Metadata_Batch)
ECMsEC <- as.fixed(Econtrols$Metadata_ECMs)
TimeEC <- as.fixed(Econtrols$Metadata_Time)
# RowEC <- as.random(Econtrols$Metadata_Row)
# ColumnEC <- as.random(Econtrols$Metadata_Column)

Econtlogmod <- aov(log(Mean_Object_Counts)~StiffnessEC + StiffnessEC%in%BatchEC +ECMsEC + TimeEC + TimeEC%in%BatchEC + StiffnessEC*TimeEC%in%BatchEC + 
                    # ECMsC%in%RowC + ECMsC%in%ColumnC+ 
                    StiffnessEC*ECMsEC + StiffnessEC*TimeEC + TimeEC*ECMsEC + 
                    StiffnessEC*ECMsEC*TimeEC, data = Econtrols)
summary(Econtlogmod)

EStiffnessNH <- as.fixed(ds_e_phenonocont$Metadata_Stiffness)
EBatchNH <- as.random(ds_e_phenonocont$Metadata_Batch)
ESceneNH <- as.random(ds_e_phenonocont$Metadata_Scene)
ERatioNH <- as.fixed(ds_e_phenonocont$Metadata_Ratio)
EECMsNH <- as.fixed(ds_e_phenonocont$Metadata_ECMs)
ETimeNH <- as.fixed(ds_e_phenonocont$Metadata_Time)
# RowNH <- as.random(ds_s_phenonocont$Metadata_Row)
# ColumnNH <- as.random(ds_s_phenonocont$Metadata_Column)

eno100logmod <- aov(log(Mean_Object_Counts)~EStiffnessNH + EStiffnessNH%in%EBatchNH + ERatioNH +EECMsNH + ETimeNH + ETimeNH%in%EBatchNH + ETimeNH*EStiffnessNH%in%EBatchNH +
                     # RatioNH%in%RowNH + RatioNH%in%ColumnNH + ECMsNH%in%RowNH + ECMsNH%in%ColumnNH+ SceneNH +
                     EStiffnessNH*ERatioNH + EStiffnessNH*EECMsNH + EStiffnessNH*ETimeNH + ETimeNH*EECMsNH + ERatioNH*ETimeNH + ERatioNH*EECMsNH + 
                     EStiffnessNH*ERatioNH*ETimeNH + EStiffnessNH*ERatioNH*EECMsNH + EStiffnessNH*EECMsNH*ETimeNH + ERatioNH*ETimeNH*EECMsNH +
                     EStiffnessNH*ERatioNH*ETimeNH*EECMsNH, data = ds_e_phenonocont)
summary(eno100logmod, digits=6)


# Graphs ------------------------------------------------------------------
Day0nohundo <- subset(ds_e_phenonocont, Metadata_Time == "Day 0", select=c(Metadata_Batch, Metadata_Stiffness, Mean_Object_Counts, Metadata_Ratio, Metadata_ECMs))
Day2nohundo <- subset(ds_e_phenonocont, Metadata_Time == "Day 2", select=c(Metadata_Batch, Metadata_Stiffness, Mean_Object_Counts, Metadata_Ratio, Metadata_ECMs))

Day0nohundo$Metadata_Ratio <- as.factor(Day0nohundo$Metadata_Ratio)
Day0nohundo$Metadata_Ratio <- factor(Day0nohundo$Metadata_Ratio, levels=c("25:75", "\"50:50\"", "\"75:25\""))

Day2nohundo$Metadata_Ratio <- as.factor(Day2nohundo$Metadata_Ratio)
Day2nohundo$Metadata_Ratio <- factor(Day2nohundo$Metadata_Ratio, levels=c("25:75", "\"50:50\"", "\"75:25\""))


Day0nohundoheatmap <-ggplot(data = Day0nohundo,
                        aes(x = Metadata_Ratio,
                            y = Metadata_ECMs,
                            fill = Mean_Object_Counts)) +
  facet_wrap(vars(Metadata_Stiffness)) +
  theme(axis.title.x = element_text(size= 30),
        axis.text.x = element_text(angle = 0, size= 20),
        axis.title.y = element_text(size= 30),
        axis.text.y = element_text(size= 15),
        strip.text =  element_text(size = 30),
        plot.title = element_text(size = 30))+
  geom_tile() +
  xlab("Ratio") +
  ylab("ECMs")+
  labs(title = "Number of Cells on Islands Day 0")
ggsave(filename = "plots/Day0HeatMapNo100.png",
       width = 15, height=10)

Day2nohundoheatmap <- ggplot(data = Day2noHOutliersRemoved,
                            aes(x = Metadata_Ratio,
                                y = Metadata_ECMs,
                                fill = Mean_Object_Counts)) +
  facet_wrap(vars(Metadata_Stiffness)) +
  theme(axis.title.x = element_text(size= 30),
        axis.text.x = element_text(angle = 0, size= 20),
        axis.title.y = element_text(size= 30),
        axis.text.y = element_text(size= 15),
        strip.text =  element_text(size = 30),
        plot.title = element_text(size = 30))+
  geom_tile() +
  xlab("Ratio") +
  ylab("ECMs")+
  labs(title = "Number of Cells on Islands Day 2")
ggsave(filename = "plots/Day2HeatMapNo100.png",
       width = 15, height= 10)

LogDay0nohundoheatmap <-ggplot(data = Day0nohundo,
                            aes(x = Metadata_Ratio,
                                y = Metadata_ECMs,
                                fill = log(Mean_Object_Counts))) +
  facet_wrap(vars(Metadata_Stiffness)) +
  theme(axis.title.x = element_text(size= 30),
        axis.text.x = element_text(angle = 0, size= 20),
        axis.title.y = element_text(size= 30),
        axis.text.y = element_text(size= 15),
        strip.text =  element_text(size = 30),
        plot.title = element_text(size = 30))+
  geom_tile() +
  xlab("Ratio") +
  ylab("ECMs")+
  labs(title = "Log(Number of Cells) on Islands Day 0")
ggsave(filename = "plots/LogDay0HeatMapNo100.png",
       width = 15, height=10)

LogDay2nohundoheatmap <- ggplot(data = Day2noHOutliersRemoved,
                             aes(x = Metadata_Ratio,
                                 y = Metadata_ECMs,
                                 fill = log(Mean_Object_Counts))) +
  facet_wrap(vars(Metadata_Stiffness)) +
  theme(axis.title.x = element_text(size= 30),
        axis.text.x = element_text(angle = 0, size= 20),
        axis.title.y = element_text(size= 30),
        axis.text.y = element_text(size= 15),
        strip.text =  element_text(size = 30),
        plot.title = element_text(size = 30))+
  geom_tile() +
  xlab("Ratio") +
  ylab("ECMs")+
  labs(title = "Log(Number of Cells) on Islands Day 2")
ggsave(filename = "plots/LogDay2HeatMapNo100.png",
       width = 15, height= 10)

Day2noHOutliersRemoved <-subset(Day2nohundo, Mean_Object_Counts < 700,select=c(Metadata_Batch, Metadata_Stiffness, Mean_Object_Counts, Metadata_Ratio, Metadata_ECMs))

ControlHeatMap <- ggplot(data= Econtrols,
                         aes(x= Metadata_Time,
                             y = Metadata_ECMs,
                             fill = Mean_Object_Counts))+
  facet_wrap(vars(Metadata_Stiffness)) +
  theme(axis.title.x = element_text(size= 30),
        axis.text.x = element_text(angle = 0, size= 20),
        axis.title.y = element_text(size= 30),
        axis.text.y = element_text(size= 15),
        strip.text =  element_text(size = 30),
        plot.title = element_text(size = 30))+
  geom_tile() +
  xlab("Time in Culture") +
  ylab("ECMs")+
  labs(title = "Number of Cells on Control Islands")
ggsave(filename = "plots/100HeatMap.png",
       width = 15, height= 10)

LogControlHeatMap <- ggplot(data= Econtrols,
                         aes(x= Metadata_Time,
                             y = Metadata_ECMs,
                             fill = log(Mean_Object_Counts)))+
  facet_wrap(vars(Metadata_Stiffness)) +
  theme(axis.title.x = element_text(size= 30),
        axis.text.x = element_text(angle = 0, size= 20),
        axis.title.y = element_text(size= 30),
        axis.text.y = element_text(size= 15),
        strip.text =  element_text(size = 30),
        plot.title = element_text(size = 30))+
  geom_tile() +
  xlab("Time in Culture") +
  ylab("ECMs")+
  labs(title = "Log(Number of Cells) on Control Islands")
ggsave(filename = "plots/log100HeatMap.png",
       width = 15, height= 10)

ggplot(data = ds_e_pheno,
       aes(x = Metadata_Stiffness,
           y = log(Mean_Object_Counts))) +
  theme(axis.title.x = element_text(size= 15),
        axis.text.x = element_text(angle = 0, size= 20),
        axis.title.y = element_text(size= 15),
        axis.text.y = element_text(size= 15),
        strip.text =  element_text(size = 15),
        plot.title = element_text(size = 15))+
geom_boxplot()

ggplot(data = ds_e_pheno,
       aes(x = Metadata_Stiffness,
           y = log(Mean_Object_Counts))) +
  facet_wrap(vars(Metadata_Time)) +
  theme(axis.title.x = element_text(size= 15),
        axis.text.x = element_text(angle = 0, size= 20),
        axis.title.y = element_text(size= 15),
        axis.text.y = element_text(size= 15),
        strip.text =  element_text(size = 15),
        plot.title = element_text(size = 15))+
  geom_boxplot()

ggplot(data = ds_e_pheno,
       aes(x = Metadata_Ratio,
           y = log(Mean_Object_Counts))) +
  theme(axis.title.x = element_text(size= 15),
        axis.text.x = element_text(angle = 0, size= 20),
        axis.title.y = element_text(size= 15),
        axis.text.y = element_text(size= 15),
        strip.text =  element_text(size = 15),
        plot.title = element_text(size = 15))+
  geom_boxplot() 
  
ggplot(data = ds_e_pheno,
       aes(x = log(Mean_Object_Counts),
           y = Metadata_ECMs)) +
  theme(axis.title.x = element_text(size= 15),
        axis.text.x = element_text(angle = 0, size= 20),
        axis.title.y = element_text(size= 15),
        axis.text.y = element_text(size= 5),
        strip.text =  element_text(size = 15),
        plot.title = element_text(size = 15))+
  geom_boxplot()

ggplot(data = ds_e_pheno,
       aes(x = log(Mean_Object_Counts),
           y = Metadata_ECMs)) +
  facet_wrap(vars(Metadata_Stiffness))+
  theme(axis.title.x = element_text(size= 15),
        axis.text.x = element_text(angle = 0, size= 20),
        axis.title.y = element_text(size= 15),
        axis.text.y = element_text(size= 5),
        strip.text =  element_text(size = 15),
        plot.title = element_text(size = 15))+
  geom_boxplot()

ggplot(data = ds_e_pheno,
       aes(x = Metadata_Ratio,
           y = log(Mean_Object_Counts))) +
  facet_wrap(vars(Metadata_Stiffness)) +
  theme(axis.title.x = element_text(size= 15),
        axis.text.x = element_text(angle = 0, size= 10),
        axis.title.y = element_text(size= 15),
        axis.text.y = element_text(size= 5),
        strip.text =  element_text(size = 15),
        plot.title = element_text(size = 15)) +
  geom_boxplot()

ggplot(data = ds_e_pheno,
       aes(x = log(Mean_Object_Counts),
           y = Metadata_ECMs)) +
  facet_wrap(vars(Metadata_Ratio))+
  theme(axis.title.x = element_text(size= 15),
        axis.text.x = element_text(angle = 0, size= 20),
        axis.title.y = element_text(size= 15),
        axis.text.y = element_text(size= 5),
        strip.text =  element_text(size = 15),
        plot.title = element_text(size = 15))+
  geom_boxplot()

ggplot(data = ds_e_pheno,
       aes(x = log(Mean_Object_Counts),
           y = Metadata_ECMs)) +
  facet_wrap(vars(Metadata_Time))+
  theme(axis.title.x = element_text(size= 15),
        axis.text.x = element_text(angle = 0, size= 20),
        axis.title.y = element_text(size= 15),
        axis.text.y = element_text(size= 5),
        strip.text =  element_text(size = 15),
        plot.title = element_text(size = 15))+
  geom_boxplot()
