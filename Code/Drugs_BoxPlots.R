suppressMessages(library(tidyr))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(ggplot2))

# Box Plot - Positive Mode ------------------------------------------------
# INPUT:
df <- fread("input/20200225_DarkWebDrug_Pos.csv", sep=",", header=TRUE)

# subset data
df <- subset(df, df$ATTRIBUTE_SampleType == "reference material")

df$ATTRIBUTE_Info_Label_Manufacturer <- factor(df$ATTRIBUTE_Info_Label_Manufacturer, 
                                               levels = c('CIPLA LTD','Combitic Global Caplet PVT. LTD.','Healing Pharma',
                                                          'fortune health care','Laborate Pharmaceuticals India LTD.','pharose remedies Ltd.',
                                                          'Naman Pharma Drugs','Signature Phytochemical Industries',
                                                          'Skohind Labs strategic alliance with Rousel Laboratories',
                                                          'uni medicolabs','Medreich limited','Centurion Laboratories PVT. LTD.'))

index_man <- fread("input/renumbering_manuf_index.txt", sep="\t", header=TRUE)

df <- merge(index_man, df, by.x="manufacturer", by.y="ATTRIBUTE_Info_Label_Manufacturer")

# PLOT the GNPS cluster Index of interest

#Octabenzone
plot_this_feature <- df$`1230_327.200779196779_2.12621626262626`

# Box
box <- ggplot(data=df, aes(x=as.factor(Cbind_Sample_ID), y=as.numeric(plot_this_feature)))+
  facet_grid(~number_man, scales = "free_x", space = "free_x")+
  geom_boxplot(aes(colour=ATTRIBUTE_API, fill=ATTRIBUTE_API), alpha=0.3, width = 0.25, outlier.alpha = 0)+
  geom_jitter(aes(fill=ATTRIBUTE_API), colour="black", position = position_jitterdodge(jitter.width = 0.5), pch=16, cex=0.7)+
  scale_color_brewer(palette = "Set1")+
  scale_fill_brewer(palette = "Set1")+
  theme_minimal()+
  theme(panel.background = element_rect(fill="grey85"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_line(colour="grey90",size=0.5, linetype="dashed"),
        axis.ticks=element_line(colour ="black",size=0.5, linetype="solid"),
        axis.text=element_text(colour="black",size=6),
        axis.text.x=element_text(angle=90, hjust=0, vjust=0.25),
        axis.line=element_line(colour="black",size=0.5, linetype="solid"),                  
        axis.title=element_text(colour="black",size=6),
        strip.text.x = element_text(colour="black",size=6),
        legend.position="none") +
  labs(x="Formulation", y="Peak Area")
print(box)
ggsave(box, filename="R_output/Drugs_Box_Pos_Octobenzone_1230_327.200779196779_2.12621626262626.pdf", width = 7.1, height = 2, units = "in", dpi=300, useDingbats=FALSE)

#2,4,7,9-Tetramethyl-5-decyne-4,7-diol	
plot_this_feature <- df$`1452_209.190190212232_4.52943333333333`

# Box
box <- ggplot(data=df, aes(x=as.factor(Cbind_Sample_ID), y=as.numeric(plot_this_feature)))+
  facet_grid(~number_man, scales = "free_x", space = "free_x")+
  geom_boxplot(aes(colour=ATTRIBUTE_API, fill=ATTRIBUTE_API), alpha=0.3, width = 0.25, outlier.alpha = 0)+
  geom_jitter(aes(fill=ATTRIBUTE_API), colour="black", position = position_jitterdodge(jitter.width = 0.5), pch=16, cex=0.7)+
  scale_color_brewer(palette = "Set1")+
  scale_fill_brewer(palette = "Set1")+
  theme_minimal()+
  theme(panel.background = element_rect(fill="grey85"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_line(colour="grey90",size=0.5, linetype="dashed"),
        axis.ticks=element_line(colour ="black",size=0.5, linetype="solid"),
        axis.text=element_text(colour="black",size=6),
        axis.text.x=element_text(angle=90, hjust=0, vjust=0.25),
        axis.line=element_line(colour="black",size=0.5, linetype="solid"),                  
        axis.title=element_text(colour="black",size=6),
        strip.text.x = element_text(colour="black",size=6),
        legend.position="none") +
  labs(x="Formulation", y="Peak Area")
print(box)

ggsave(box, filename="R_output/Drugs_Box_2,4,7,9-Tetramethyl-5-decyne-4,7-diol_1452_209.190190212232_4.52943333333333.pdf", width = 7.1, height = 2, units = "in", dpi=300, useDingbats=FALSE)


# Box Plot - Negative Mode ------------------------------------------------
# INPUT:
df <- fread("input/20200225_DarkWebDrug_Neg.csv", sep=",", header=TRUE)

# subset data
df <- subset(df, df$ATTRIBUTE_SampleType == "reference material")

df$ATTRIBUTE_Info_Label_Manufacturer <- factor(df$ATTRIBUTE_Info_Label_Manufacturer, 
                                               levels = c('CIPLA LTD','Combitic Global Caplet PVT. LTD.','Healing Pharma',
                                                          'fortune health care','Laborate Pharmaceuticals India LTD.','pharose remedies Ltd.',
                                                          'Naman Pharma Drugs','Signature Phytochemical Industries',
                                                          'Skohind Labs strategic alliance with Rousel Laboratories',
                                                          'uni medicolabs','Medreich limited','Centurion Laboratories PVT. LTD.'))

index_man <- fread("input/renumbering_manuf_index.txt", sep="\t", header=TRUE)

df <- merge(index_man, df, by.x="manufacturer", by.y="ATTRIBUTE_Info_Label_Manufacturer")

# PLOT the GNPS cluster Index of interest

# putative SDS related
plot_this_feature <- df$`72_265.148674994159_4.78577658486708`

# Box
box <- ggplot(data=df, aes(x=as.factor(Cbind_Sample_ID), y=as.numeric(plot_this_feature)))+
  facet_grid(~number_man, scales = "free_x", space = "free_x")+
  geom_boxplot(aes(colour=ATTRIBUTE_API, fill=ATTRIBUTE_API), alpha=0.3, width = 0.25, outlier.alpha = 0)+
  geom_jitter(aes(fill=ATTRIBUTE_API), colour="black", position = position_jitterdodge(jitter.width = 0.5), pch=16, cex=0.7)+
  scale_color_brewer(palette = "Set1")+
  scale_fill_brewer(palette = "Set1")+
  theme_minimal()+
  theme(panel.background = element_rect(fill="grey85"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_line(colour="grey90",size=0.5, linetype="dashed"),
        axis.ticks=element_line(colour ="black",size=0.5, linetype="solid"),
        axis.text=element_text(colour="black",size=6),
        axis.text.x=element_text(angle=90, hjust=0, vjust=0.25),
        axis.line=element_line(colour="black",size=0.5, linetype="solid"),                  
        axis.title=element_text(colour="black",size=6),
        strip.text.x = element_text(colour="black",size=6),
        legend.position="none") +
  labs(x="Formulation", y="Peak Area")
print(box)
ggsave(box, filename="R_output/Drugs_Box_Neg_unknownputativeSDS_72_265.148674994159_4.78577658486708.pdf", width = 7.1, height = 2, units = "in", dpi=300, useDingbats=FALSE)

# putative SDS related
plot_this_feature <- df$`84_293.180258751828_5.34812832310838`

# Box
box <- ggplot(data=df, aes(x=as.factor(Cbind_Sample_ID), y=as.numeric(plot_this_feature)))+
  facet_grid(~number_man, scales = "free_x", space = "free_x")+
  geom_boxplot(aes(colour=ATTRIBUTE_API, fill=ATTRIBUTE_API), alpha=0.3, width = 0.25, outlier.alpha = 0)+
  geom_jitter(aes(fill=ATTRIBUTE_API), colour="black", position = position_jitterdodge(jitter.width = 0.5), pch=16, cex=0.7)+
  scale_color_brewer(palette = "Set1")+
  scale_fill_brewer(palette = "Set1")+
  theme_minimal()+
  theme(panel.background = element_rect(fill="grey85"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_line(colour="grey90",size=0.5, linetype="dashed"),
        axis.ticks=element_line(colour ="black",size=0.5, linetype="solid"),
        axis.text=element_text(colour="black",size=6),
        axis.text.x=element_text(angle=90, hjust=0, vjust=0.25),
        axis.line=element_line(colour="black",size=0.5, linetype="solid"),                  
        axis.title=element_text(colour="black",size=6),
        strip.text.x = element_text(colour="black",size=6),
        legend.position="none") +
  labs(x="Formulation", y="Peak Area")
print(box)
ggsave(box, filename="R_output/Drugs_Box_Neg_unknownputativeSDSanalog_84_293.180258751828_5.34812832310838.pdf", width = 7.1, height = 2, units = "in", dpi=300, useDingbats=FALSE)
